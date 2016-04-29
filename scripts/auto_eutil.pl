#!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::EUtilities;
use Getopt::Long;
use Bio::SeqIO;

my $email = '';
my $logging;
my $quiet = 0;
my $manualfile;

#Get options
my $signal = GetOptions(
                         'log:s'    => \$logging,
                         'quiet'    => \$quiet,
                         'manual=s' => \$manualfile,
                         'email=s'  => \$email
                       );

if ( !$email ) {
    die("Must provide e-mail address with -email. REQUIRED by NCBI for submission to their servers.\n"
       );
}

if ( defined($logging) ) {
    if ( $logging eq '' ) {
        $logging = "eutil.log";
    }
}

my $infile = shift;

open INFILE, "$infile"
  or die "$infile is unavailable. Check your path and try again!\n";
my @ids = <INFILE>;
close INFILE;
chomp(@ids);
for (@ids) {

#This removes the revision number from the accession, if present. Required for mapping GIs later
    s/\.[\d]+//;
}
my %idmatch;    # master record accession for WGS; accession for nt
my %accn;       # key = gi num, value = acc num
my ( %assem, @assem );    # key = gi num, value = assembly id
my ( %taxon, @taxon );    # key = gi num, value = taxon id
my (%names);              # key = taxon id, value = scientific name
my ( %gis, @gis );        # key = acc num, value = gi num
my (%guess);              # key = gi num, value = best guess at species name
my @cc = culture();       # Culture collections, used for guessing names later
my %cc;                   #Map culture collections to keys in hash for checking
my %titles;
my %host;
my %country;
my %source;
my %gb_name;
my %strain;
my %culture;

for (@cc) {
    $cc{$_} = 1;
}
my @link_chunks;          # value = reference to each chunk of ids
my $chunk_size = 500;
my $iter       = 0;
my $chunk_num  = 0;
my $history;    # a cookie with WebEnv and query_key for later retrieval;
                #Bio::Tools::EUtilities::History object
my $i = 0;
my %counter;    # key = name; value = counts of retrieved records
my %manual;     # print these IDs to manual file
my %memory;     # Keep track of genome names to identify conflicts


#Must chunk the data before submitting to EUtils

foreach my $id (@ids) {

    #Convert wgs records to base record to find assembly ids
    my $master = getMaster($id);
    if ( !exists $idmatch{$master} ) {
        $idmatch{$master} = $id;
    }
}

foreach my $master ( sort keys %idmatch ) {
    if ( $iter < $chunk_size ) {
        push( @{ $link_chunks[$chunk_num] }, $master );
        $iter++;
    } else {
        $iter = 0;
        $chunk_num++;
        push( @{ $link_chunks[$chunk_num] }, $master );
    }
}

#Initialize EUtil factory
#Efetch to retrieve GI numbers from accessions

my $factory = Bio::DB::EUtilities->new(
                                        -eutil => 'efetch',
                                        -db    => 'nucleotide',
                                        -email => $email,
                                      );

my $j = 1;

logger( "Submitting " . ( scalar(@link_chunks) ) . " chunks to efetch.\n" );

$counter{'assem'}    = 0;    # Number of assemblies retrieved
$counter{'taxon'}    = 0;    # Number of names retrieved
$counter{'nuc'}      = 0;    # Number of GIs retrieved
$counter{'organism'} = 0;    # Number of organism names retrieved

foreach my $chunk (@link_chunks) {

    #Submit each chunk, retrieve GI numbers as glob
    $factory->set_parameters( -id      => $chunk,
                              -rettype => 'gi' );
    my @gi_chunk;

    if ( @gi_chunk = split( m{\n}, $factory->get_Response->content() ) ) {
        logger("Genome id chunk $j successfully fetched\n");
    }

    $counter{'nuc'} += scalar(@gi_chunk);    #count number of GIs retrieved
    push( @gis, \@gi_chunk )
      ;    #Save references to each chunk in @gis, will submit to epost next
    $j++;
}

logger( $counter{'nuc'} . " GIs collected using efetch.\n" );

#Set parameters for epost, history to use the GIs on both esummary and elink
#Correspondence to relate GI to result

$factory->set_parameters( -eutil => 'epost', );

logger( "Submitting " . scalar(@gis) . " chunks to epost.\n" );

$j = 1;

foreach my $gi_chunk (@gis) {

    $factory->reset_parameters(
                                -eutil          => 'epost',
                                -db             => 'nucleotide',
                                -email          => $email,
                                -id             => $gi_chunk,
                                -keep_histories => 1,
                              );

    if ( $history = $factory->next_History ) {
        logger(   "GI chunk $j submitted to epost successfully, with "
                . scalar(@$gi_chunk)
                . " IDs submitted\n" );
    } else {
        logger(
             "Problem submitting GI chunk $j to epost. Unrecoverable error.\n");
        die;
    }

    #Pass history objects to eutil factory to retrieve gis posted using epost
    #Use esummary to retrieve assembly and taxon ids from nucleotide summary
    #Also link GI to accession using the summary
    $factory->set_parameters( -eutil   => 'esummary',
                              -history => $history );

    logger("Retrieving taxon ids\n");

    while ( my $ds = $factory->next_DocSum ) {

        my $id = $ds->get_id;    # $id = gi number

        while ( my $item = $ds->next_Item('flattened') ) {
            my $name = $item->get_name;
            if ( $name eq 'Caption' ) {

                #Accession numbers stored in Caption of esummary xml file
                my $acc = $item->get_content;
                if ( !$acc ) {
                    logger("Unable to find accession number for GI $id\n");
                    next;
                }
                $gis{$acc} = $id;
                $accn{$id} = $acc;
            }
            if ( $name eq 'TaxId' ) {

                #Taxon ids stored in TaxId of esummary xml file
                my $taxid = $item->get_content;
                if ( !$taxid ) {
                    logger("Unable to find taxon id for GI $id\n");
                    next;
                }
                $taxon{$id} = $taxid;
                $counter{'taxon'}++;
            }
            if ( $name eq 'Title' ) {

     #Title includes the name of the sequence, and potentially the genus/species
                my $title = $item->get_content;
                my ( $genus, $species, @split ) = split( ' ', $title );
                my $sciname = join( " ", $genus, $species );
                if ( $sciname =~ /,/ ) {
                    $sciname =~ s/,//;
                    $guess{$id} = $sciname;
                    next;
                }
                for ( my $k = 0 ; $k < scalar(@split) ; $k++ ) {
                    my $value = $split[$k];
                    if ( exists( $cc{$value} ) ) {
                        $sciname .= " $split[$k] $split[$k+1]";
                        $sciname =~ s/,//;
                        last;
                    }
                    if (    $value =~ /contig/i
                         || $value =~ /chromosome/i
                         || $value =~ /scaffold/i
                         || $value =~ /linear/i
                         || $value =~ /circular/i
                         || $value =~ /plasmid/i
                         || $value =~ /complete/i
                         || $value =~ /DNA/
                         || $value =~ /^main$/
                         || $value =~ /whole/i
                         || $value =~ /draft/i )
                    {
                        #Add more values as they are found.
                        last;
                    }
                    if ( $value =~ /,/ ) {
                        $value =~ s/,//;
                        $sciname .= " $value";
                        last;
                    }
                    $sciname .= " $value";
                }
                $guess{$id} = $sciname;
            }
        }
    }
    my $retry    = 0;
    my $retstart = 0;
    my $retmax   = 500;

    my $count = @$gi_chunk;

    logger(
        "Retrieving Host, Isolation Source, and Country data, if available!\n");

    $factory->set_parameters(
                              -eutil     => 'efetch',
                              -history   => $history,
                              -db        => 'nuccore',
                              -seq_start => 1,
                              -seq_stop  => 2,
                              -rettype   => 'gbwithparts'
                            );

  RETRIEVE_SEQS:
    while ( $retstart < $count ) {
        $factory->set_parameters( -retmax   => $retmax,
                                  -retstart => $retstart );
        my @gbfile;
        eval { @gbfile = split( m{\n}, $factory->get_Response->content ); };
        if ($@) {
            die "Server error: $@.  Try again later" if $retry == 5;
            print STDERR "Server error, redo #$retry\n";
            $retry++ && redo RETRIEVE_SEQS;
        }

        my $gbstring = join( "\n", @gbfile ) . "\n";

        open my $stringfh, "<", \$gbstring
          or die "Unable to open string for reading: $!";

        my $seqio = Bio::SeqIO->new( -fh     => $stringfh,
                                     -format => 'genbank' );

        while ( my $gb_seq = $seqio->next_seq() ) {
            my $gb_acc   = $gb_seq->accession;
            my $organism = 'NA';
            my $strain   = 'NA';
            my $host     = 'NA';
            my $country  = 'NA';
            my $note     = 'NA';
            my $source   = 'NA';
            my $culture  = 'NA';
            my $isolate  = 'NA';
            foreach my $feat ( $gb_seq->get_SeqFeatures ) {

                if ( $feat->primary_tag eq "source" ) {
                    if ( $feat->has_tag('organism') ) {
                        my @data = $feat->get_tag_values('organism');
                        $organism = $data[0];
                    }
                    if ( $feat->has_tag('strain') ) {
                        my @data = $feat->get_tag_values('strain');
                        $strain = $data[0];
                    }
                    if ( $feat->has_tag('host') ) {
                        my @data = $feat->get_tag_values('host');
                        $host = $data[0];
                    }
                    if ( $feat->has_tag('country') ) {
                        my @data = $feat->get_tag_values('country');
                        $country = $data[0];
                    }
                    if ( $feat->has_tag('note') ) {
                        my @data = $feat->get_tag_values('note');
                        $note = $data[0];
                    }
                    if ( $feat->has_tag('isolation_source') ) {
                        my @data = $feat->get_tag_values('isolation_source');
                        $source = $data[0];
                    }
                    if ( $feat->has_tag('culture_collection') ) {
                        my @data = $feat->get_tag_values('culture_collection');
                        $culture = $data[0];
                        $culture =~ s/:/ /;
                    }
                    if ( $feat->has_tag('isolate') ) {
                        my @data = $feat->get_tag_values('isolate');
                        $isolate = $data[0];
                    }
                }
            }

            #print STDERR "DEBUG => ".join("\t",'id',$gb_gi,'acc',$gb_acc,
            #                              'organism',$organism,'strain',
            #                              $strain,'host',$host,'country',
            #                              $country,'note',$note) . "\n";
            my $gb_master = getMaster($gb_acc);

            my $gb_gi = $gis{$gb_master};

            if ( !$gb_gi ) {
                print STDERR
                  "ERROR: Unable to find gi num for $gb_acc (MASTER=$gb_master)\n";
            }

            my $gb_name = 'NA';

            if ( $organism =~ /$strain/ ) {
                $gb_name = $organism;
            } else {
                if ( $strain ne 'NA' ) {
                    $gb_name = join( " ", $organism, $strain );
                } elsif ( $culture ne 'NA' ) { 
                    $gb_name = join( " ", $organism, $culture );
                } elsif ( $isolate ne 'NA' ) {
                    $gb_name = join( " ", $organism, $isolate );
                    $strain = $isolate;
                }

            }

            $host{$gb_gi}    = $host;
            $country{$gb_gi} = $country;
            $source{$gb_gi}  = $source;
            $gb_name{$gb_gi} = $gb_name;
            $strain{$gb_gi}  = $strain;
            $culture{$gb_gi} = $culture;

        }
        $retstart += $retmax;
    }


    #Set eutil factory to retrieve links from nucleotide GIs to assembly IDs
    $factory->reset_parameters(
                                -eutil          => 'elink',
                                -history        => $history,
                                -dbfrom         => 'nuccore',
                                -db             => 'assembly',
                                -cmd            => 'acheck',
                                -email          => $email,
                                -keep_histories => 1
                              );

    logger("Getting links to assembly database from nucleotide GI numbers\n");
    my %assemlinks = ();
    my %retry      = ();

    # iterate through the LinkSet objects
    while ( my $ds = $factory->next_LinkSet ) {
        my $submit = ( $ds->get_submitted_ids )[0];
        my @items  = $ds->get_link_names;
        my $match  = grep( /assembly/, @items );
        if ($match) {
            $assemlinks{$submit} = 1;
            next;
        } else {

#If no match is found, get accession of original record (not master) and retry with that.
            my $acc = $accn{$submit};
            if ( $acc eq $idmatch{$acc} ) {
                $assem{$submit} = 'NULL';
            } else {
                $retry{ $idmatch{$acc} } = 1;
            }
        }
    }

    if (%retry) {
        my @retry = sort keys %retry;
        my $retry_factory =
          Bio::DB::EUtilities->new(
                                    -eutil   => 'efetch',
                                    -db      => 'nucleotide',
                                    -email   => $email,
                                    -rettype => 'gi',
                                    -id      => \@retry
                                  );

        #        $factory->reset_parameters(
        #                                    -eutil   => 'efetch',
        #                                    -id      => \@retry,
        #                                    -rettype => 'gi',
        #                                    -db      => 'nucleotide',
        #                                    -email   => $email
        #                                  );
        my @ret_chunk;

        if ( @ret_chunk =
             split( m{\n}, $retry_factory->get_Response->content() ) )
        {
            logger(   "Retrying assembly links for "
                    . scalar(@retry)
                    . " accessions...\n" );
        }

        $retry_factory->set_parameters( -eutil => 'epost', );

        $retry_factory->set_parameters( -id             => \@ret_chunk,
                                        -keep_histories => 1, );
        my $rethistory;
        if ( $rethistory = $retry_factory->next_History ) {
            logger(
                   "Retry for GI chunk $j submitted to epost successfully with "
                     . scalar(@ret_chunk)
                     . " GIs submitted\n" );
        }

        $retry_factory->set_parameters( -eutil   => 'esummary',
                                        -history => $rethistory );

        while ( my $ds = $retry_factory->next_DocSum ) {

            my $id = $ds->get_id;    # $id = gi number

            while ( my $item = $ds->next_Item('flattened') ) {
                my $name = $item->get_name;
                if ( $name eq 'Caption' ) {

                    #Accession numbers stored in Caption of esummary xml file
                    my $acc = $item->get_content;
                    if ( !$acc ) {
                        logger("Unable to find accession number for GI $id\n");
                        next;
                    }
                    $gis{$acc} = $id;
                    $accn{$id} = $acc;
                }
            }

        }
        $retry_factory->set_parameters(
                                        -eutil  => 'elink',
                                        -dbfrom => 'nuccore',
                                        -db     => 'assembly',
                                        -cmd    => 'acheck',
                                        -id     => $rethistory
                                      );

        while ( my $ds = $retry_factory->next_LinkSet ) {
            my $submit = ( $ds->get_submitted_ids )[0];
            my @items  = $ds->get_link_names;
            my $match  = grep( /assembly/, @items );
            if ($match) {
                $assemlinks{$submit} = 1;
                next;
            } else {

                $assem{$submit} = 'NULL';
            }
        }

    }

    my @assemlinks = keys %assemlinks;

    my $assem_fac =
      Bio::DB::EUtilities->new(
                                -eutil          => 'elink',
                                -db             => 'assembly',
                                -dbfrom         => 'nuccore',
                                -correspondence => 1,
                                -email          => $email,
                                -cmd            => 'neighbor',
                                -id             => \@assemlinks
                              );

    logger("Saving assembly links found for GI chunk $j\n");

    while ( my $ds = $assem_fac->next_LinkSet ) {
        my $submit  = ( $ds->get_submitted_ids )[0];
        my @results = $ds->get_ids;
        if ( scalar(@results) > 1 ) {
            logger(
                "GI Number: $submit has multiple assembly results.  Check to make sure the proper name is in the final keyfile\n"
            );
        }
        my $result = pop(@results);

        my $acc    = $accn{$submit};
        my $master = getMaster($acc);

        $assem{$master} = $result;
        $counter{'assem'}++;
    }
    $j++;
}

logger("$counter{'assem'} assembly links found.\n");
logger("$counter{'taxon'} taxonomy links found.\n");

@taxon = values(%taxon);    #Store all taxids in @taxon

#logger("Taxon ids:\n".join("\n",@taxa));

#Chunk taxon IDs to post to NCBI

my @tax_chunks;
$iter      = 0;
$chunk_num = 0;

foreach my $taxon (@taxon) {
    if ( $iter < $chunk_size ) {
        push( @{ $tax_chunks[$chunk_num] }, $taxon );
        $iter++;
    } else {
        $iter = 0;
        $chunk_num++;
        push( @{ $tax_chunks[$chunk_num] }, $taxon );
    }
}

logger( "Submitting " . ( scalar(@tax_chunks) ) . " chunks to epost.\n" );

$j = 1;

#Set EUtil factory to submit using epost to taxonomy database

foreach my $chunk (@tax_chunks) {

    #Submit each chunk of taxids
    $factory->set_parameters(
                              -eutil          => 'epost',
                              -db             => 'taxonomy',
                              -id             => $chunk,
                              -keep_histories => 1,
                              -use_history    => 'n'
                            );

    if ( $history = $factory->next_History ) {
        logger("Taxon id chunk $j posted to epost successfully\n");
    } else {
        logger("Unable to fetch next taxon history. Unrecoverable error.\n");
        die;
    }

    #Use epost histories on esummary to get taxonomy summaries for each taxid
    $factory->set_parameters(
                              -eutil       => 'esummary',
                              -use_history => 'y',
                              -history     => $history
                            );

    logger("Getting taxonomic names.\n");
    while ( local $_ = $factory->next_DocSum ) {
        my $taxonid  = ( $_->get_contents_by_name('TaxId') )[0];
        my @sci_name = $_->get_contents_by_name('ScientificName');
        my $name     = pop(@sci_name);
        $names{$taxonid} =
          $name;    # $taxonid is each taxid; value is each scientific name

    }

    $j++;
}

#Map accessions to taxonomy names
foreach my $id (@ids) {
    my $master = getMaster($id);
    my (
         $gi,   $taxid,   $name,   $assemid, $guess,
         $host, $country, $source, $gb_name, $strain,
         $culture
       )
      = getData($master);

    if ( $gi eq 'NULL' ) {
        next;
    }
    if ( $guess ne 'NULL' ) {
        if ($manualfile) {
            my @split  = split( " ", $guess );
            my @split2 = split( " ", $name );
            my $file   = $infile;
            $file =~ s/.fas.accn.tmp/.key/;
            my $words = 2;
            if ( $guess =~ /candidatus/i || $name =~ /candidatus/i ) {
                $words = 3;
            }
            if ( scalar(@split) == $words && scalar(@split2) == $words ) {
                logger(
                    "Check $manualfile for $id, needs manual examination to get appropriate species name\n"
                );
                open my $fh, ">>", "$manualfile"
                  or die "Unable to open file for manual checking! : $!";
                print $fh join( "\t", $file, $id, $guess, $guess ) . "\n";
                close $fh;
            }
        }
    }

    print join( "\t",
                $id,     $assemid, $taxid, $name,    $guess,  $gi,
                $master, $gb_name, $host,  $country, $source, $strain,
                $culture
              )
      . "\n";

}

logger("Done.\n");

sub getMaster {
    my $id = shift;
    if ( $id =~ /[A-Z]{4}_?[0-9]{8}/ ) {
        $id =~ s/[0-9]/0/g;
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

sub culture {
    my @cc = (
               'LMG',   'ATCC',  'DSM',   'JCM',   'BCCM', 'BCRC',
               'BKM',   'VKM',   'CDC',   'NCDC',  'CECT', 'CFBP',
               'CIP',   'DSMZ',  'KCTC',  'KCCM',  'NBRC', 'IFO',
               'NCIB',  'NCIMB', 'NCPPB', 'ECCO',  'CCUG', 'MCC',
               'HAMBI', 'NCCB',  'JGI',   'CCBAU', 'NRRL', 'KFB'
             );
    return (@cc);
}

sub getData {
    my $master = shift;
    my $gi;         #gi number, from %gis, $acc is key
    my $taxid;      #taxonomy id, from %taxon, $gi is key
    my $name;       #scientific name, from %names, $taxid is key
    my $assemid;    #assembly id, $gi is key
    my $guess;      #best guess at name, $gi is key
    my $host;
    my $country;
    my $source;
    my $gb_name;
    my $strain;
    my $culture;

    $gi = $gis{$master};
    if ( !$gi ) {
        logger(
            "Unable to find GI number for accession $master. This accession will not be included in the final analysis.\n"
        );
        $gi = 'NULL';
    }
    $assemid = $assem{$master};
    if ( !$assemid ) {
        $assemid = 'NULL';
    }
    $taxid = $taxon{$gi};
    if ($taxid) {
        $name = $names{$taxid};
        if ( !$name ) {
            $name = $master;
        }
    } else {
        $taxid = 'NULL';
        $name  = $master;
    }
    $guess = $guess{$gi};
    if ( !$guess ) {
        $guess = 'NULL';
    }
    $host = $host{$gi};
    if ( !$host || $host eq 'NA' ) {
        $host = 'NULL';
    }
    $country = $country{$gi};
    if ( !$country || $country eq 'NA' ) {
        $country = 'NULL';
    }
    $source = $source{$gi};
    if ( !$source || $source eq 'NA' ) {
        $source = 'NULL';
    }
    $gb_name = $gb_name{$gi};
    if ( !$gb_name || $gb_name eq 'NA' ) {
        $gb_name = 'NULL';
    }
    $strain = $strain{$gi};
    if ( !$strain || $strain eq 'NA' ) {
        $strain = 'NULL';
    }
    $culture = $culture{$gi};
    if ( !$culture || $culture eq 'NA' ) {
        $culture = 'NULL';
    }

    return (
             $gi,   $taxid,   $name,   $assemid, $guess,
             $host, $country, $source, $gb_name, $strain,
             $culture
           );
}

sub manPrint {

}
