#!/usr/bin/env perl
#######################################################################
#
# COPYRIGHT NOTICE
#
# autoANI - A program to automate Average Nucleotide Identity calcs
#
# Copyright (C) 2015
# Edward Davis
# Jeff Chang
#
# This file is a part of autoANI.
#
# autoANI is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# autoANI is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more detail.
#
# You should have received a copy of the GNU General Public License
# along with autoANI.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################
use warnings;
use strict;
use Bio::SeqIO;
use File::Spec;
use Getopt::Long;
use Pod::Usage;
use Parallel::ForkManager;
use Time::HiRes qw{sleep};
my $stime = 0.33;

#Set Defaults here
my $size = 1020;
my $pid_cutoff = 30;
my $coverage = 70;

#Set to complete path for the blast bins if not found in PATH
my $blastbindir = '';
#Set to scripts dir if you move autoANI.pl to somewhere in your PATH
my $scriptdir   = '';

#Set to your e-mail address if you will be the primary user
my $email = '';
if (defined $ENV{'EMAIL'}) {
    $email ||= $ENV{'EMAIL'};
}

my $incommand = join( " ", $0, @ARGV );

my $cwd = File::Spec->curdir();
$cwd = File::Spec->rel2abs($cwd);
my $blastout = "$cwd/blast";
my @dbs;
my @chunked;

my $removeNs = 0;    #set to 0 to keep Ns in sequence.
my $help     = 0;
my $man      = 0;
my %headers;
my $finish = 0;
my $quiet  = 0;
my $log    = 1;
my %idmatch;
my $threads = 1;
my $retkeys;
my $outfile;
my %matched;
my $blasthitsfile;
my $sge = 0; # sge flag only works on Oregon State University CGRB infrastructure
my $blast_v_check = qr/2.2.31|2.[3-9]+.[\d]|3.[\d]+.[\d]+/;
my $prep;
my $oldani;
my $queue = 'bpp'; #Only relevant for Oregon State University CGRB infrastructure

my ( $svol, $sdir, $sfile ) = File::Spec->splitpath($0);
$sdir .= "scripts/";

my $program     = $blastbindir . 'blastn';
my $makeblastdb = $blastbindir . 'makeblastdb';
my $elinkpath   = $scriptdir . 'auto_edirect.pl';

if ( !$scriptdir ) {
    $elinkpath   = $sdir . 'auto_edirect.pl';
}

my $version = '1.3.0';
my $date = 'November 15, 2016';
my $vhelp;

my $signal = GetOptions(
                         'size=i'      => \$size,
                         'pid=i'       => \$pid_cutoff,
                         'coverage=f'  => \$coverage,
                         'help'        => \$help,
                         'man'         => \$man,
                         'finish'      => \$finish,
                         'quiet'       => \$quiet,
                         'log!'        => \$log,
                         'threads=i'   => \$threads,
                         'email=s'     => \$email,
                         'keys'        => \$retkeys,
                         'outfile=s'   => \$outfile,
                         'blasthits=s' => \$blasthitsfile,
                         'sge'         => \$sge,
                         'prep'        => \$prep,
                         'oldani'      => \$oldani,
                         'version'     => \$vhelp,
                         'queue=s'       => \$queue
                       );
#sge flag only works on Oregon State University CGRB infrastructure.

if ($help) {
    pod2usage( -verbose => 1,
               -output  => ">&STDERR" );
}

if ($man) {
    pod2usage( -verbose => 2,
               -output  => ">&STDERR" );
}

if ($vhelp) {
    pod2usage( -verbose => 0,
               -output  => ">&STDERR",
               -msg     => "$0 version $version.\nReleased on $date.\n");
}

die("Unknown option passed.  Check parameters and try again.\n") if !$signal;

if ( $sge ) {
    my $hostname = `hostname`;
    chomp($hostname);
    if ( $hostname !~ /cgrb\.oregonstate\.local/ ) {
        logger("You do not appear to be on the Oregon State University infrastructure.\n");
        logger("Setting \$sge flag to 0\n");
        $sge = 0;
    } else {
        my @submit_hosts = `qconf -ss`;
        chomp(@submit_hosts);
        my $success = 0;
        foreach my $submit_host (@submit_hosts) {
            if ( $hostname =~ /$submit_host/ ) {
                logger("$hostname is SGE submit host.\n");
                $success = 1;
                last;
            }
        }
        if ( $success == 0 ) {
            logger(
                "$hostname is not submit host.\n"
                . "Retry submission from a submit host below "
                . "(or run the script locally without the -sge flag):\n"
                . join( "\n", @submit_hosts )
                . "\n" 
            );
            exit(-1);
        }
    }   
}

my @infiles;
foreach (@ARGV) {
    my $infile = File::Spec->rel2abs($_);
    push( @infiles, $infile );
}

my %infiles = ();

%infiles = map { $_ => 1 } @infiles if (@infiles);

if ( -s "./queries/queries.txt" ) {
    open my $queries, "<", "./queries/queries.txt"
      or die
      "Unable to open queries.txt (from ./queries/queries.txt). Check permissions and try again.\n";
    my @queries = <$queries>;
    chomp(@queries);
    foreach my $query (@queries) {
        $infiles{$query} = 1;
    }
    close $queries;
}

die("No input specified!") if ( !%infiles );

if ( !$email ) {
    die("Must provide a valid email address to continue!  Set the email within the file so that you do not have to provide it on the command line each time."
       );
}

if ( $pid_cutoff > 70 ) {
    logger(
        "Percent identity cutoff set above 70%. Must be between 0 and 70. Setting to 70.\n"
    );
    $pid_cutoff = 70;
} elsif ( $pid_cutoff < 0 ) {
    logger(
        "Percent identity cutoff set below 0%. Must be between 0 and 70. Setting to 0.\n"
    );
    $pid_cutoff = 0;
}

if ( $coverage > 100 ) {
    logger(
        "Coverage cutoff set to above 100. Must be between 0 and 100. Setting to 100.\n"
    );
} elsif ( $coverage < 0 ) {
    logger(
        "Coverage cutoff set to below 0. Must be between 0 and 100. Setting to 0.\n"
    );
}

my $time = localtime();
logger("Command as submitted at $time:\n$incommand\n");

my @blast_version = `$program -version`;

if ( grep( /$blast_v_check/, @blast_version ) ) {
    logger("Found BLAST version 2.2.31 or greater\n");
} else {
    logger(
        "BLAST version 2.2.31 or greater is REQUIRED!  Make sure the full path is provided or include it in your PATH!\n"
    );
    exit(-1);
}


my @makeblastdb_check = `$makeblastdb -version`;
logger(join('',@makeblastdb_check));
if ( $? != 0 ) {
    logger("Problem generating blast databases!\n");
    logger("Make sure makeblastdb is in your PATH or you have added the path to the top of $0\n");
    exit(-1);
}

#my $makeblastdb_check = system($makeblastdb,"-version");
#
#if ($makeblastdb_check) {
#    logger("Unable to find makeblastdb.\n");
#    logger("Check your system settings and try again.\n");
#    exit(-1);
#} else {
#    logger("Found makeblastdb.\n");
#}

#Setup blastdbs
my @folders = ( "db", "fasta", "blast", "queries" );

foreach my $folder (@folders) {
    if ( !-d $folder ) {
        my $error = system("mkdir ./$folder");
        die "Unable to generate $folder folder! Check writing permissions. : $!"
          if $error != 0;
    }
}

logger("Saving query paths to ./queries/queries.txt\n");

open my $queries, ">", "./queries/queries.txt"
  or die
  "Unable to open queries.txt for writing. Check permissions and try again!\n";
@infiles = sort keys %infiles;
print $queries join( "\n", @infiles ) . "\n";
close $queries;

foreach my $infile1 ( sort keys %infiles ) {
    my ( $vol1, $dir1, $file1 ) = File::Spec->splitpath($infile1);
    foreach my $infile2 ( sort keys %infiles ) {
        next if $infile1 eq $infile2;
        my ( $vol2, $dir2, $file2 ) = File::Spec->splitpath($infile2);
        if ( $file1 eq $file2 ) {
            logger("Found indentically named file in two locations:\n");
            logger("$infile1\n");
            logger("$infile2\n\n");
            logger(
                "Please check these paths and rename one of the files (if appropriate), or remove the duplicate file.\n"
            );
            logger(
                "You must remove the offending file from ./queries/queries.txt to continue.\n"
            );
            exit(-1);
        }
    }
    if ( !-s $infile1 ) {
        logger("Unable to find input file $infile1.\n");
        logger(
            "Check the path name and correct it in ./queries/queries.txt before continuing.\n"
        );
        exit(-1);
    }
}

$time = localtime();

logger("Setting up blast databases and splitting FASTA files at $time\n");

open my $keyfile, ">>", "ani.keys" or die "Unable to open keyfile  : $!";

foreach my $infile (@infiles) {
    my ( $volume, $dir, $file ) = File::Spec->splitpath($infile);
    push( @chunked, "$cwd/fasta/$file.$size" );
    my $in = Bio::SeqIO->new( -file   => "$infile",
                              -format => 'fasta' );
    my $skip = 0;
    if ( -e "fasta/$file.$size" ) {
        if ( -s "fasta/$file.$size" ) {
            $skip = 1;
            logger("Already found $file.$size. Skipping\n");
        } else {
            `rm -f fasta/$file.$size`;
        }
    }
    if ( $skip == 0 ) {
        open my $temp, ">>", "ani.accn.tmp"
          or die "Unable to open ani.accn.tmp : $!";
        logger("Generating chunked file $file.$size from $file.\n");
        open my $fastaout, ">", "fasta/$file.$size"
          or die "$cwd/fasta/$file.$size is unavailble : $!";
        my $count = 0;
        eval {
            while ( my $seq = $in->next_seq() ) {
                my $name  = $seq->id;
                my $i     = 0;
                my $chunk = 0;
                my $desc  = $seq->description;
                my $local = 0;

                get_keys( $name, $desc, $temp );

                if ( $removeNs == 1 ) {
                    my $sequence = $seq->seq;
                    $sequence =~ s/N//g;
                    $seq->seq($sequence);
                }
                my $length = $seq->length;
                while ( $i < $length ) {
                    my $min = $i + 1;
                    my $max = 0;
                    if ( ( $i + $size ) > $length ) {
                        last;
#                        $max = $length;
                    } else {
                        $max = $i + $size;
                    }
                    print $fastaout ">${name}_${chunk} $desc\n";
                    print $fastaout $seq->subseq( $min, $max ) . "\n";
                    $i += $size;
                    $chunk++;
                }
                $count += $chunk;
            }
        };

        if ($@) {
            `rm -f ./fasta/$file.$size`;
            logger("There was a problem parsing $file. Check to ensure it is a FASTA formatted file and try again.\n");
            logger("*****If $file is not a FASTA formatted file, then it needs to be removed from ./queries/queries.txt before continuing.*****\n");
            logger("Error Message:\n" . $@);
            exit(-1);
        }
        close $fastaout;

        if ( $count == 0 ) {
            logger("Unable to chunk file $file.\n");
            logger("Check to ensure this file is a FASTA formatted file.\n");
            exit(-1);
        }

        if ( ! -s "./fasta/$file.$size" ) {
            logger("There was a problem generating the chunked file ./fasta/$file.$size.\n");
            logger("Check to ensure your input file was in FASTA format and try again.\n");
            exit(-1);
        }
    }

    if (    ( -e "db/$file.nhr" )
         && ( -e "db/$file.nin" )
         && ( -e "db/$file.nsq" ) )
    {
        logger("Blast database already present for db/$file. Skipping...\n");
        push( @dbs, "$cwd/db/$file" );
        next;
    } else {
        my $output = `$makeblastdb -in $infile -dbtype nucl -out db/$file`;
        logger( $output . "\n" );
        push( @dbs, "$cwd/db/$file" );
        if ( $? != 0 ) {
            logger("Problem generating blast databases!\n");
            exit(-1);
        }
    }
}

$time = localtime();

logger("Done at $time.\n");

if ($retkeys) {
    logger("Refinding keys from headers.\n");
    `rm -f ani.accn.tmp`;
    `rm -f ani.keys`;
    open my $temp, ">>", "ani.accn.tmp"
      or die "Unable to open ani.accn.tmp : $!";
    foreach my $infile (@infiles) {
        open my $headers, "-|",
          "grep '>' $infile  | cut -c 2- "
          or die "Unable to get FASTA headers : $!\n";
        my @headers = <$headers>;
        s/_[0-9]+$// for @headers;
        my %headers = map { $_ => 1 } @headers;

        foreach my $header ( sort keys %headers ) {
            my ( $name, $desc ) = split( " ", $header, 2 );
            get_keys( $name, $desc, $temp );
        }
    }
    close $temp;
}

close $keyfile;

if ( -s "ani.accn.tmp" ) {
    $time = localtime();
    logger("Retrieving taxonomic data from NCBI at $time\n");
    my $command = join( " ",
                        $elinkpath, "--log",           "ani.eutil.log",
                        "--email",
                        $email,     "ani.accn.tmp",    ">>",
                        "ani.keys" );
    my $error = system("$command");
    if ( $error != 0 ) {
        logger(
             "Unable to retrieve taxnomic data! Resubmit command and try again!"
        );
        die("\n");
    }
    `rm -f ani.accn.tmp`;
    my @sorted = `cat ani.keys | sort | uniq`;
    `rm -f ani.keys`;
    open my $sortfh, ">", "ani.keys" or die "Unable to open ani.keys! $!\n";
    print $sortfh @sorted;
    $time = localtime();
    logger("Done at $time\n");
}

my $calcs = @infiles;

$calcs = $calcs * ($calcs - 1);

if ( $prep ) {
    logger("$calcs BLAST searches required before completion.\n");
    logger("-prep flag given. Quitting before BLAST.\n");
    exit();
}

if ( $finish == 0 ) {
    logger("Starting $calcs BLAST searches.\n");

    %headers = ();

    my $outfmt = q{6 qseqid sseqid pident length evalue nident};

    my $pm = Parallel::ForkManager->new( $threads - 1 );

    $time = localtime();
    logger("BLAST searches started at $time\n");

    foreach my $querypath (@chunked) {
        my ( $vol, $dir, $query ) = File::Spec->splitpath($querypath);
        foreach my $dbpath (@dbs) {
            $pm->start and next;
            my ( $vol1, $dir1, $subject ) = File::Spec->splitpath($dbpath);
            my $test = $query;
            $test =~ s/\.$size//;
            if ( $test ne $subject ) {
                my $output = "./blast/";
                $output .= $query . "_vs_" . $subject . ".tab";
                if ( -s $output ) {
                    logger("BLAST output $output already found!\n");
                } else {
                    if ( -e $output ) {
                        `rm -f $output`;
                    }
                    #@command = join( " ",
                    #                 $blastscript, $coverage, $size,
                    #                 $pid_cutoff, $output, @command );

                    if ( $sge == 0 ) {
                        my @command =       ($program,                 
                                            "-db",              qq{$dbpath},    
                                            "-query",           qq{$querypath}, 
                                            "-outfmt",          qq{$outfmt},                  
                                            "-evalue",          "0.001",   
                                            "-num_threads",     1,         
                                            "-max_target_seqs", "1",      
                                            "-max_hsps",        "1",      
                                            "-task",            "blastn", 
                                            "-xdrop_gap",       150,      
                                            "-penalty",         "-1",     
                                            "-dust",            "no",
                                            "-reward",          1,        
                                            "-gapopen",         5,        
                                            "-gapextend",       2,
                                            "-out",             $output);
                        logger("Running BLAST command:\n@command\n");
                        my $check = system(@command);
                        if ( $check != 0 || ! -s $output ) {
                            `rm -f $output`;
                            logger("BLAST search between $query and $subject failed.\n");
                            logger("Check your input files and try again.\n");
                            exit(-1)
                        }

                    } else {
                        #Only works on Oregon State University CGRB Infrastructure
                        my @command =       ($program,                 
                                            "-db",              qq{$dbpath},    
                                            "-query",           qq{$querypath}, 
                                            "-outfmt",          qq{\'$outfmt\'},                  
                                            "-evalue",          "0.001",   
                                            "-num_threads",     1,         
                                            "-max_target_seqs", "1",      
                                            "-max_hsps",        "1",      
                                            "-task",            "blastn", 
                                            "-xdrop_gap",       150,      
                                            "-penalty",         "-1",     
                                            "-dust",            "no",
                                            "-reward",          1,        
                                            "-gapopen",         5,        
                                            "-gapextend",       2,
                                            "-out",             $output);
                        logger("Running BLAST command through SGE:\n@command\n");
                        @command = join( " ",
                                         "SGE_Batch",                "-r",
                                         "sge.${query}_vs_$subject", "-q",
                                         $queue,                      "-c",
                                         qq{"@command"},             "-Q" );
                        my @runoutput = `@command`;
                        sleep($stime);
                    }
                }

            }
            $pm->finish;
        }
    }
    $pm->wait_all_children;

    $time = localtime();
    logger("BLAST searches completed at $time\n");

    logger("To complete analysis, use the -finish flag.\n");
    exit(0);
}

#Concatenate data to analyze
logger("Calculating results from BLAST output...\n");

#`rm -f ./blast/all.out` if -e "./blast/all.out";
#my $command = "cat ./blast/*.tab > ./blast/all.out";
#$error = system($command);
#if ( $error != 0 ) {
#    `rm -f ./blast/all.out` if -e "./blast/all.out";
#    die "Problem concatenating output files!\n";
#}

logger("Loading names from ani.keys\n");

if ( !-s "ani.keys" ) {
    keys_warn();
}

open my $keys, "-|", "sort ani.keys | uniq"
  or die "ani.keys is unavailable : $!\n";

#Current format is Accession,AssemblyID,TaxID,SciName,GuessName,GI,Master,GenBankName,Host,Country,Source,Strain,CultureCollection
#Current format is Accession,AssemblyID,TaxID,SciName,GI,Master,GenBankName,Country,Source,Strain,CultureCollection,Year
#                  $match        0        1      2    3    4         5         6       7      8          9           10
while (<$keys>) {
    my $line = $_;
    chomp($line);
    my ( $accn, @values ) = split( "\t", $line );
    my $match = $values[0];
    if ( $match eq 'NULL' ) {
        $match = $values[4];
    }
    $idmatch{$accn} = $match;

    if ( !exists( $headers{$match} ) ) {
        
        #Choose appropriate header name
        my $sciname = $values[2];
        my $gbname  = $values[5];
        my $header  = $sciname;
        my @test = split( " ", $sciname );  #Test for quality of name from taxid
        my $test = @test;
        my @test2 = split( " ", $gbname );
        my $test2 = @test2;
        my $candidatus = 0;
        my $subsp = 0;
        if ($sciname =~ /candidatus/i ) {
            $candidatus = 1;
        }
        if ($sciname =~ /pv\.|bv\.|subsp\./i ) {
            $subsp = 2;
        }
        my $value = 2 + $candidatus + $subsp;
        if ( $test == $value ) {
            if ( $gbname ne 'NULL' ) {
                if ( $sciname ne $gbname ) {
                    $header = $gbname;
                }
            }
        }

        $headers{$match} = $header;
    }

}
close $keys;

logger("Printing name information to ani_rename.log\n");

open my $rename, ">", "ani_rename.log"
  or die "Unable to open ani_rename.log : $!\n";

foreach my $match ( sort keys %headers ) {
    print $rename join( "\t", $headers{$match}, $match ) . "\n";
}

close $rename;

if ( !%idmatch ) {
    keys_warn();
}

$time = localtime();

logger("Extracting BLAST output started at $time\n");

my %results;
my %tmpresults;
my %genomes;
my $tmpfile = 'ani.out.tmp';
my %seenfh;
my %filenames;

if ( -s $tmpfile ) {
    logger("Temp file ($tmpfile) found. Extracting previous results now.\n");
    open my $tmpfh, "<", "$tmpfile" or die "Unable to open $tmpfile : $!\n";
    while(<$tmpfh>) {
        my $line = $_;
        chomp($line);
        my ($query, $subject, $ani, $ani_new, $hits, $file) = split("\t",$line);
        if ($file) {
            if (-s "./blast/$file") {
                $seenfh{$file} = 1;
            }
        }
        $tmpresults{$query}{$subject}{'ani'} = $ani;
        $tmpresults{$query}{$subject}{'ani_new'} = $ani_new;
        $tmpresults{$query}{$subject}{'hits'} = $hits;
        $genomes{$query} = 1;
    }
}

my @blastout = `find ./blast/ -type f -exec basename {} \\;`;

chomp(@blastout);

my $blastfiles = @blastout;

if ($blastfiles < 2) {
    logger("Unable to find BLAST output in ./blast/\n");
    logger("Check input files and ./queries/queries.txt to ensure proper files are being examined.\n");
    exit(-1);
}

$time = localtime();

logger("Found $blastfiles BLAST output files to parse. Starting at $time.\n");

foreach my $boutfile (sort @blastout) {
    next if $seenfh{$boutfile};

    open my $blast, "<", "./blast/$boutfile"
      or die "Unable to open blast output file $boutfile : $!\n";

#Analyze data

    while (<$blast>) {
        my $line = $_;
        chomp($line);
        next if ( $line =~ /#/ );
        my @data      = split( "\t", $line );
        my $query     = $data[0];
        my $subject   = $data[1];
        my $percentid = $data[2];
        my $length    = $data[3];
        my $nident    = $data[5];
        my $new_pid   = ( $nident / $size ) * 100;
        my ( $query_contig, $subject_contig );
        my ( $query_genome, $subject_genome );
        $query =~ s/_[0-9]+(?=$)//;

        my ( $qacc, $qmatched ) = get_accn($query);
        my ( $sacc, $smatched ) = get_accn($subject);

        my @accessions = ( $qacc, $sacc );

        foreach my $accession (@accessions) {
            if ( !defined( $idmatch{$accession} ) ) {
                logger(
                    "Unable to find assembly/identifier for accession $accession. Check ani.keys to make sure this value is found.\n"
                );
                exit(-1);
            }
        }

        $query_genome   = $idmatch{$qacc};
        $subject_genome = $idmatch{$sacc};

        $genomes{$query_genome} = 1;
        $filenames{$query_genome}{$subject_genome} = $boutfile;

        if ( $query ne $subject ) {
            next
              if (    $length < ( ( $coverage / 100 ) * $size )
                   || $new_pid < $pid_cutoff );
            push( @{ $results{$query_genome}{$subject_genome} }, [ $percentid, $new_pid ] );
        }
    }

    close $blast;
}
#
#    open my $tmpfile, "<", "ani.tmp"
#      or die "Unable to find temp file ani.tmp. Cannot calculate results!\n";
#
#    while (<$tmpfile>) {
#        my $line = $_;
#        chomp($line);
#        my ( $qacc, $sacc, $ani, $hits ) = split( "\t", $line );
#        my @accessions = ( $qacc, $sacc );
#
#        foreach my $accession (@accessions) {
#            if ( !defined( $idmatch{$accession} ) ) {
#                logger(
#                    "Unable to find assembly/identifier for accession $accession. Check ani.keys to make sure this value is found.\n"
#                );
#                exit(-1);
#            }
#        }
#
#        my $query_genome   = $idmatch{$qacc};
#        my $subject_genome = $idmatch{$sacc};
#
#        $genomes{$query_genome} = 1;
#
#        $tmpresults{$query_genome}{$subject_genome}{'hits'} = $hits;
#        $tmpresults{$query_genome}{$subject_genome}{'ani'}  = $ani;
#
#    }
#}
#
$time = localtime();

logger("Finished parsing BLAST output at $time\n\n");

my $outfh;
my $hitsfh;

if ($outfile) {
    open $outfh, ">", $outfile
      or die "Unable to open output file $outfile : $!";
} else {
    $outfh = \*STDOUT;
}

my @genomes = sort( keys(%genomes) );
my @print;
foreach my $genome (@genomes) {
    push( @print, $headers{$genome} );
}
print $outfh join( "\t", "", @print ) . "\n";

if ($blasthitsfile) {
    open $hitsfh, ">", $blasthitsfile
      or die "Unable to open output hits file $blasthitsfile: $!";
    print $hitsfh join( "\t", "", @print ) . "\n";
}

logger("Attempting to print output for " . @genomes . " genomes.\n");

my $countoutfh = 0;

foreach my $genome (@genomes) {
    print $outfh $headers{$genome} . "\t";
    if ($blasthitsfile) {
        print $hitsfh $headers{$genome} . "\t";
    }
    foreach my $genome2 (@genomes) {
        if ( $genome eq $genome2 ) {
            print $outfh "----------\t";
            print $hitsfh "----------\t" if $blasthitsfile;
            next;
        }

        my $ani = undef;
        my $ani_new = undef;
        my $n   = 0;

        if ( defined( $results{$genome}{$genome2} ) ) {
            my $sum = 0;
            my $sum_new = 0;
            my $n   = scalar( @{ $results{$genome}{$genome2} } );
            foreach my $result ( @{ $results{$genome}{$genome2} } ) {
                $sum += $result->[0];
                $sum_new += $result->[1];
            }
            $ani = sprintf( "%.3f", $sum / $n );
            $ani_new = sprintf( "%.3f", $sum_new / $n );
            open my $tmpfh, ">>", "ani.out.tmp" or die "Unable to open ani.out.tmp file : $!";
            print $tmpfh join("\t",$genome,$genome2,$ani,$ani_new,$n,$filenames{$genome}{$genome2}) . "\n";
            close $tmpfh;
        } elsif ( defined( $tmpresults{$genome}{$genome2} ) ) {
            $ani = sprintf( "%.3f", $tmpresults{$genome}{$genome2}{'ani'} );
            $ani_new = sprintf( "%.3f", $tmpresults{$genome}{$genome2}{'ani_new'} );
            $n = $tmpresults{$genome}{$genome2}{'hits'};
        }

        if ($ani) {
            $countoutfh++;
        } else {
            logger("Unable to find ANI value for $genome to $genome2 comparison.\n");
            $ani = 'NA';
        }
        my $printani = $ani_new;
        if ($oldani) {
            $printani = $ani;
        }
        print $outfh $printani . "\t";
        print $hitsfh $n . "\t" if $blasthitsfile;
    }
    print $outfh "\n";
    print $hitsfh "\n" if $blasthitsfile;
}
$time = localtime();
logger("Printed output for $countoutfh comparisons.\n");
logger("Finished at $time.\n");

sub logger {
    my $message = shift;
    print STDERR $message unless $quiet == 1;
    if ( $log == 1 ) {
        open my $fh, ">>", "ani.log" or die "ani.log is unavailable : $!";
        print $fh $message;
        close $fh;
    }
}

sub get_accn {
    my $name = shift;
    my @split = split( '\|', $name );
    my @accn_types = (
                qr/\A[A-Z][0-9]{5}(\.[0-9])?\Z/,
                qr/\A[A-Z]{2}_?[0-9]{6}(\.[0-9])?\Z/,
                qr/\A([A-Z]{2}_)?[A-Z]{4}[0-9]{8}([0-9]{1,2})?(\.[0-9])?\Z/,
                qr/\A[A-Z]{2}_[A-Z]{2}[0-9]{6}(\.[0-9])?\Z/
    );
    my $match;
    my $matched = 0;
    foreach my $data (@split) {
        foreach my $type (@accn_types) {
            if ( $data =~ $type ) {
                $match = $data;
                $matched = 1;
            }
        }
    }
    if (!$match) {
        if ($split[1]) {
            $match = $split[1];
        } else {
            $match = $name;
        }
    }
    return ($match,$matched);
}

#sub get_accn {
#    my $name = shift;
#    my @split = split( '\|', $name );
#    if ( $name !~ /gnl/ ) {
#        my @accn_types = (
#                    qr/\A[A-Z][0-9]{5}(\.[0-9])?\Z/,
#                    qr/\A[A-Z]{2}_?[0-9]{6}(\.[0-9])?\Z/,
#                    qr/\A([A-Z]{2}_)?[A-Z]{4}[0-9]{8}([0-9]{1,2})?(\.[0-9])?\Z/,
#                    qr/\A[A-Z]{2}_[A-Z]{2}[0-9]{6}(\.[0-9])?\Z/
#        );
#        foreach my $data (@split) {
#            foreach my $type (@accn_types) {
#                if ( $data =~ $type ) {
#                    return $data;
#                }
#            }
#
#        }
#    } else {
#        return $split[1];
#    }
#}

sub keys_warn {
    logger(
        "Unable to find names in ani.keys.  Check to make sure ani.keys was properly generated.  If ani.keys is empty, delete the files in the ./fasta folder and resubmit the command.\n"
    );
    exit();
}

sub get_keys {
    my $name = shift;
    my $desc = shift;
    my $temp = shift;
    my ($accession, $matched) = get_accn($name);
    if ($matched) {
        if ( !$accession ) {
            logger(
                "Unable to find accession for $name. Check proper format and try again.\n"
            );
            die;
        }
        print $temp $accession . "\n";
    } else {
        if ( !$accession ) {
            logger(
                "Unable to find accession for local genome $name. Check proper format and try again.\n"
            );
            die;
        }
        $desc = $1 if $desc =~ /\[([^\]+])\]$/;
        $desc =~ s/\[//;
        $desc =~ s/\]//;
        my $title = '';
        if ($desc) {
            $title = $desc;
        } else {
            logger("Unsure if file containing $name is formatted correctly. Check your input to be sure you submitted the correct file.\n");
            $title = $accession;
        }
        print $keyfile join( "\t",
                             $accession, 'NULL', 'NULL', $title,
                             $title,     'NULL', $accession, $title,
                             'NULL',     'NULL', 'NULL',     'NULL',
                             'NULL')
          . "\n";
    }
}
__END__

=head1 NAME

autoANI.pl - Calculate pairwise ANI (Average Nucleotide Identity) values between any number of sequenced genomes

=head1 SYNOPSIS

autoANI.pl input1.fasta input2.fasta ... input[n].fasta

=head1 OPTIONS

Defaults shown in square brackets.  Possible values shown in parentheses.

=over 8

=item B<-help|h>

Print a brief help message and exits.

=item B<-man>

Print a verbose help message and exits.

=item B<-quiet>

Turns off progress messages. Use noquiet to turn messages back on.

=item B<-version>

Print version number and exit.

=item B<-email> (email@univ.edu) - B<REQUIRED>

Enter your email.  NCBI requires this information to continue.

=item B<-log|nolog> [logging on]

Using -nolog turns off logging. Logfile is in the format of ani.log.

=item B<-threads> [1]

*Recommended* - Using multiple threads will significantly speed up the BLAST searches.

=item B<-prep>

Do the prep and quit.

=item B<-finish>

Finish the job after BLAST searches.

=item B<-oldani>

Calculate the ANI value using the 'old' algorithm, where there was no adjustment for alignment length.

=item B<-keys>

Re-generate keyfile from NCBI/local FASTA files.  Generally unnecessary unless ani.keys has been altered/removed.

=item B<-outfile> 

Optional.  Path to output file. By default, output is directed to STDOUT.

=item B<-size> [1020]

Genome chunk size.

=item B<-coverage> [70] (0-100)

Percentage of query coverage cutoff.

=item B<-pid> [30] (0-70)

Percent identity cutoff for BLAST search results.

=back

=head1 DESCRIPTION

This program expects to take in one multifasta/fasta file per genome (incomplete genomes are fine) and chunks the genome into the specified sizes (1020 default).  The chunked sequences are then used as queries in BLAST searches against the other full genomes provided.  After removing BLAST results below the cutoff (70% coverage and 30% ID), the average identity of the rest is calculated.  Output is a tab-delimited file with all of the pairwise comparisons.

This script handles any type of nucleotide file from the NCBI nucleotide database (WGS included) as long as the accession number is provided in the header line of the FASTA file.  See examples below.

>U00096.3
>BA000007.2
>gi|254160123|ref|NC_012967.1|

The script also handles FASTA files generated using fna_rename.pl (with the following format)

>gnl|K-12_substr_MDS42 complete genome [Escherichia coli K-12 substr MDS42]

The important parts are the >gnl|(Unique identifier) and the [Genus species strain] at the end.

=cut
