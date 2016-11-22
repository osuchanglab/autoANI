#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my $edirpath = ''; #Set /path/to/edirect/ if edirect tools are not in your path
my $elink = $edirpath . 'elink';
my $epost = $edirpath . 'epost';
my $efetch = $edirpath . 'efetch';
my $xtract = $edirpath . 'xtract';

my $email = '';
my $logging;
my $quiet = 0;

#Get options
my $signal = GetOptions(
                         'log:s'    => \$logging,
                         'quiet'    => \$quiet,
                         'email=s'  => \$email
                       );

if ( defined($logging) ) {
    if ( $logging eq '' ) {
        $logging = "eutil.log";
    }
}

checkEmail($email);

my $infile = shift;
my @ids;

open INFILE, "$infile"
  or die "$infile is unavailable. Check your path and try again!\n";
while(<INFILE>) {
    my $line = $_;
    chomp($line);
    if ($line) {
        push(@ids,$line);
    }
}
close INFILE;

my %idmatch;    # master record accession for WGS; accession for nt
my @idmatch;    # master record list for submit to NCBI
my %accn;       # key = gi num, value = acc num
my ( %assem, @assem );    # key = gi num, value = assembly id
my ( %taxon, @taxon );    # key = gi num, value = taxon id
my (%names);              # key = taxon id, value = scientific name
my ( %gis, @gis );        # key = acc num, value = gi num
my %titles;
my %country;
my %source;
my %gb_name;
my %strain;
my %culture;
my $history;
my %counter;
my %year;

my %memory;     # Keep track of genome names to identify conflicts

foreach my $id (@ids) {

    #Convert wgs records to base record to find assembly ids
    my $master = getMaster($id);
    if ( !exists $idmatch{$master} ) {
        $idmatch{$master} = $id;
    }
}

@idmatch = sort keys %idmatch;

#Initialize EUtil factory
#Efetch to retrieve GI numbers from accessions

my $j = 1;

logger( "Submitting " . ( scalar(@idmatch) ) . " ids to epost.\n" );


my $join = join(',',@idmatch);

if ( -e 'history' ) {
    `rm -f history`;
}

my $command = join(" ", $epost, '-format', 'acc', '-db', 'nuccore', '-id', "$join", '>', 'history');

system("$command") == 0 or die "Unable to run epost command : $!";

if ( $? != 0 ) {
    logger("Problem posting accession numbers to NCBI. See output message below and check your connection and try again.\n");
    exit(4);
}

$history = `cat history`;

$counter{'epost'} = getCount($history);

logger("Successfully posted $counter{'epost'} accessions.\n");

#This section will become extraneous once accession.versions are acceptable.

@gis = `cat history | $efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion Gi`;

chomp(@gis);

logger( scalar(@gis) . " GIs returned from efetch.\n" );

if (scalar(@gis) == scalar(@idmatch)){
    #print STDERR "Returned equal number of GIs as submitted accessions!\n";
    logger("Matching Accessions to GI numbers.\n");
    foreach my $line (@gis) {
        my ($acc, $gi) = split("\t",$line);
        $accn{$gi} = $acc;
        $gis{$acc} = $gi;
    }
} else {
    logger("GIs and accessions do not match!\n");
    logger("Send your test dataset to davised.dev\@gmail.com for further troubleshooting.\n");
    logger("Script will continue, but output is probably missing data.\n");
}

#This section will remain mostly unchanged

logger("Getting genbank metadata for organism name, country, year, etc.\n");
my @command = ( "cat history |",
                "$efetch -format gbc -mode xml -seq_start 1 -seq_stop 2 |",
                "$xtract -insd source organism strain isolation_source db_xref country culture_collection collection_date");

$command = join(" ", @command);
my @output = `$command`;

logger("Parsing genbank metadata.\n");

foreach my $item (@output) {
    chomp($item);
    my ($acc,$organism,$strain,$source,$xref,$country,$cc,$date) = split("\t",$item);
    my $gi = $gis{$acc};
    if (! $gi) {
        logger("Unable to find GI for accession $acc. Skipping.\n");
        next;
    }
    $cc =~ s/:/ / if $cc =~ /:/;
    my $taxid;
    if ($xref =~ /taxon/) {
        my $junk;
        ($junk, $taxid) = split(":",$xref);
    } else {
        $taxid = '-';
    }
    my $gb_name;
    if ( $organism =~ /$strain/ ) {
        $gb_name = $organism;
    } else {
        if ( $strain ne '-' ) {
            $gb_name = join( " ", $organism, $strain );
        } elsif ( $cc ne '-' ) { 
            $gb_name = join( " ", $organism, $cc );
        } 
    }

    $names{$gi} = $organism;
    $country{$gi} = $country;
    $source{$gi}  = $source;
    $gb_name{$gi} = $gb_name;
    $strain{$gi}  = $strain;
    $culture{$gi} = $cc;
    $taxon{$gi} = $taxid;
    $year{$gi} = $date;
}

logger("Getting links to assembly DB.\n");
#This section will use accession instead of GI numbers
@command = ("cat history |",
            "$elink -target assembly -cmd neighbor |",
            "$xtract -pattern LinkSet -element IdList/Id Link/Id");
$command = join(" ", @command);

@output = `$command`;

logger("Parsing links to assembly DB.\n");

foreach my $item (@output) {
    chomp($item);
    my ($gi, $assemid) = split("\t",$item);
    $assem{$gi} = $assemid;
}

logger("Mapping and printing all data.\n");

#Map accessions to taxonomy names
foreach my $id (@ids) {
    my $master = getMaster($id);
    my (
         $gi,   $taxid,   $name,   $assemid, 
         $country, $source, $gb_name, $strain,
         $culture, $year
       )
      = getData($master);

    if ( $gi eq 'NULL' ) {
        next;
    }

    print join( "\t",
                $id,     $assemid, $taxid, $name,    $gi,
                $master, $gb_name, $country, $source, $strain,
                $culture, $year
              )
      . "\n";

}

logger("Done.\n");
`rm -f history`;

sub getMaster {
    my $id = shift;
    for (my $i = 1; $i < 20; $i++) {
        my $ipad = sprintf("%02d", $i);
        my $regex = qr/[A-Z]{4}_?$ipad/;
#        print $regex;
        if ( $id =~ $regex ) {
#            print "HERE\n";
            $id =~ s/(?<!\.)[0-9]/0/g;
            $id =~ s/\.[\d]+/.$i/;
            last;
        }
    }
    return $id;
}

sub logger {
    my $message = shift;
    print STDERR $message if $quiet == 0;

    if ($logging) {
        open my $log, ">>", "$logging" or die "Unable to open log file : $!";
        print $log $message;
        close $log;
    }
}

sub getData {
    my $master = shift;
    my $gi;         #gi number, from %gis, $acc is key
    my $taxid;      #taxonomy id, from %taxon, $gi is key
    my $name;       #scientific name, from %names, $gi is key
    my $assemid;    #assembly id, $gi is key
    my $country;
    my $source;
    my $gb_name;
    my $strain;
    my $culture;
    my $year;

    $gi = $gis{$master};
    if ( !$gi ) {
        logger(
            "Unable to find GI number for accession $master. This accession will not be included in the final analysis.\n"
        );
        $gi = 'NULL';
    }
    $assemid = $assem{$gi};
    if ( !$assemid ) {
        $assemid = 'NULL';
    }
    $taxid = $taxon{$gi};
    if (!$taxid) {
        $taxid = 'NULL';
    }
    $name = $names{$gi};
    if (!$name) {
        $name = $master;
    }
    $country = $country{$gi};
    if ( !$country || $country eq '-' ) {
        $country = 'NULL';
    }
    $source = $source{$gi};
    if ( !$source || $source eq '-' ) {
        $source = 'NULL';
    }
    $gb_name = $gb_name{$gi};
    if ( !$gb_name || $gb_name eq '-' ) {
        $gb_name = 'NULL';
    }
    $strain = $strain{$gi};
    if ( !$strain || $strain eq '-' ) {
        $strain = 'NULL';
    }
    $culture = $culture{$gi};
    if ( !$culture || $culture eq '-' ) {
        $culture = 'NULL';
    }
    $year = $year{$gi};
    if ( !$year || $year eq '-' ) {
        $year = 'NULL';
    }

    return (
             $gi,   $taxid,   $name,   $assemid,  $country, 
             $source, $gb_name, $strain, $culture, $year
           );
}

sub checkEmail {
    my $em = shift;
    if ( ! defined($ENV{'EMAIL'}) ) {
        if ( $em ) {
            if ( $em !~ /\A[\w.]+\@[\w.]+\Z/ ) {
                print STDERR "Email \'$em\' does not appear to be a valid address.\n";
                print STDERR "Check your input and try again.\n";
                exit(3);
            }
            $ENV{'EMAIL'} = $em;
        } else {
            print STDERR 'Please provide an e-mail to continue! You must set your $EMAIL environment variable for edirect to work!\n';
            exit(2);
        }
    }
}

sub getCount {
    my $epdata = shift;
    my $ret = $1 if $epdata =~ m#<Count>([\d]+)</Count>#;
    return $ret;
}

