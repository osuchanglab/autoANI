#!/usr/bin/env perl
use warnings;
use strict;
use File::Fetch;
use Getopt::Long;

my $help = 0;
my $type = ''; # Options are mlsa or ani
my $fmt  = 'abbr';
my $dl = 0;
my $debug = 0;
my $rep = 0;

if (scalar(@ARGV) == 0) {
    printHelp();
}

#Get options
my $signal = GetOptions(
                         'type=s'   => \$type,
                         'f|format=s' => \$fmt,
                         'd|download' => \$dl,
                         'h|help'   => \$help,
                         'debug'    => \$debug,
                         'r|rep'    => \$rep
                       );

my $term = shift;

my $edirectpath = ''; # Set to /path/to/edirect/ if not in $PATH
my $esearch = $edirectpath . 'esearch';
my $efetch  = $edirectpath . 'efetch';
my $xtract  = $edirectpath . 'xtract';
my $efilter = $edirectpath . 'efilter';
my $elink   = $edirectpath . 'elink';

if ($help || !$term) {
    printHelp();
}

my $email = ''; #set to email address 
if (defined $ENV{'EMAIL'}) {
    $email ||= $ENV{'EMAIL'};
}

if (!$email) {
    die "Set email parameter on line 5 of this script, or set your EMAIL environment variable to continue.";
}

if ( $fmt !~ /abbr|full|strain/ ) {
    print STDERR "Format: $fmt is invalid.\n";
    print STDERR "Possible values for output format are abbr, full, or strain.\n";
    printHelp();
    exit(-1);
}

if ($type !~ /mlsa|ani/i) {
    print STDERR "Type: '$type' is invalid or missing.\n";
    print STDERR "Valid types are mlsa or ani.\n";
    print STDERR "Please check your settings and try again.\n";
    exit(-1);
}

my $query = qq{$term};

if ( $term !~ /ORGANISM|ORGN/ ) {
    $query = qq{$query\[ORGANISM\]};
} 

#Can add other things to query here

$query = qq{'$query'};

my $command = "$esearch -db genome -query $query | $elink -batch -target assembly";

if ( $rep ) {
    $command .= " | $efilter -query 'representative\[Properties\]'";
}

my $assembly = `$command`;
my $count = getCount($assembly);
print STDERR "Found $count genomes to download.\n";
print STDERR "Expect ".($count*5)."MB to ".($count*7)."MB of data for download.\n";
if (!$dl) {
    if ($count > 0) {
        print STDERR "If you want to download the $count genomes, please resubmit with the -download flag.\n";
    } else {
        print STDERR "Your search term '$term' did not return any results. Please check that it is a valid NCBI taxonomy term and re-submit.\n";
    }
    exit(0);
} else {
    print STDERR "Starting genome downloads...\n";
}
#esearch -db genome -query 'Rhizobiaceae[ORGANISM]' | elink -target assembly | efetch -format docsum > assem_test.xml
#esearch -db genome -query 'Rhizobiaceae[ORGANISM]' | elink -target assembly | elink -target nuccore -name assembly_nuccore_insdc | efetch -format fasta > ! rhizobiaceae.fasta
#xtract -pattern DocumentSummary -element Genbank,AssemblyName,FtpPath_GenBank,Organism,Sub_value
open(my $historyfh, '>', 'assemlinks') or die "Unable to open assemlinks file : $!";
print $historyfh $assembly;
close($historyfh);

my $unknown = 1; #Counter for organisms with no strain name

if ( $type =~ /ani/i ) {
    my @command = ("cat assemlinks |",
                   "$efetch -format docsum |",
                   "$xtract -pattern DocumentSummary",
                   "-block DocumentSummary -element",
                   "Genbank",
                   "AssemblyName",
                   "FtpPath_GenBank",
                   "SpeciesName",
                   "-block DocumentSummary -element",
                   "Sub_value",
                   "-block DocumentSummary -unless",
                   q{Sub_value -lbl '\-'});
    my $command = join(" ",@command);
    my @docsums = `$command`;
    if ( $? != 0 ) {
        die("Unable to download document summaries : $!");
    }
    if ($debug == 1) {
        open(my $dfh, ">", 'docsums.txt') or die "Unable to open docsums.txt : $!";
        print $dfh @docsums;
        close($dfh);
    }
    my $i = 1;
    foreach my $line (@docsums) {
        chomp($line);
        my @data = split("\t", $line);
        my @name = split(" ", $data[3]);
        my $genus = $name[0];
        my $species = $name[1];
        my $strain = $data[4];
        if ($strain eq '-' || $data[3] =~ /$strain/) {
            if (scalar(@name) == 3) {
                $strain = $name[2];
            } elsif (scalar(@name) > 3) {
                $strain = $name[2];
                foreach my $item (@name[3..$#name]) {
                    $strain .= " $item";
                }
            } else {
                $strain = 'STRAIN'.$unknown;
                $unknown++;
            }
        }
        if ($strain) {
            $strain =~ s/ /_/g;
            $strain =~ y/\//_/;
            $strain =~ s/\(//g;
            $strain =~ s/\)//g;
            $strain =~ s/://g;
            $strain =~ s/;//g;
        }
        my $outname = '';
        if ($fmt =~ /full/ || $species =~ /sp\./) {
            $species =~ s/sp\./sp/;
            $outname = $genus.'_'.$species.'_'.$strain;
        } elsif ($fmt =~ /abbr/) {
            $outname = substr($genus, 0, 1).'_'.$species.'_'.$strain;
        } else {
            $outname = $strain;
        }
        $outname .= '.fna';
        my $gbassem = $data[0];
        my $assemname = $data[1];
        $assemname =~ s/ /_/g;
        my $ftp = $data[2];
        my $filename = join('_',$gbassem,$assemname,'genomic.fna.gz');
        my $newfile = $filename;
        $newfile =~ s/\.gz$//;
        if ($debug) {
            print STDERR "DEBUG:$newfile\t$outname\n";
            next();
        }
        if ( -s "$outname" ) {
            print STDERR "$outname already found. Skipping...\n";
            $i++;
            next;
        }
        my $uri = $ftp."/".$filename;
        my $ff = File::Fetch->new(uri => $uri);
        print STDERR "Genome $i: $filename -> $outname.\n";
        my $where = $ff->fetch() or die $ff->error;
        if ( ! -s "$where") {
            die "Unable to download file $filename. Resubmit and try again.";
        } else {
            print STDERR "Downloaded $filename. Unzipping and moving to $outname.\n";
            system("gunzip $where") == 0 or die "Unable to unzip $where using gunzip";
            system("mv $newfile $outname") == 0 or die "Unable to renamed file $newfile";
            print STDERR "Done.\n";
        }
        $i++;

        #print $outname."\n";
    }
} elsif ($type =~ /mlsa/i ) {
    # elink -target nuccore -name assembly_nuccore_insdc | efetch -format fasta > ! rhizobiaceae.fasta
    my $outname = "$term.fasta";
    if (-s "$outname") {
        print STDERR "Already found output file: $outname\n";
        print STDERR "Generate a blast database from this file using makeblastdb to use in autoMLSA.\n";
        print STDERR "Otherwise, remove or move $outname to continue.\n";
        exit(0);
    }
    my @command = ('cat assemlinks |',
                   "$elink -batch -target nuccore -name assembly_nuccore_insdc |",
                   "$efetch -format fasta > $outname");
    my $command = join(" ", @command);
    system($command) == 0 or die "Unable to download $term genome sequences!\n";
    if ( ! -s "$term.fasta" ) {
        print STDERR "Unable to find $outname file.\n";
        print STDERR "Problem downloading and/or saving the file.\n";
        print STDERR "Check your settings and permissions and try again.\n";
        exit(-1);
    } else {
        print STDERR "Downloaded genomes to $outname\n";
        print STDERR "Generate a blast database from this file using makeblastdb to include in autoMLSA pipeline.\n";
    }

}

sub getCount {
    my $c = shift;
    my $ret = $1 if $c =~ m#<Count>([\d]+)</Count>#;
    return $ret;
}

sub printHelp {
    print STDERR "Please re-submit command with this syntax:\n";
    print STDERR "$0 -type [mlsa or ani] -format [format] -rep [NCBI Taxonomy term]\n";
    print STDERR "Acceptable output formats are full, abbr (default), and strain\n";
    print STDERR "Example: $0 -type mlsa -format abbr Pseudomonas\n";
    exit();
}
