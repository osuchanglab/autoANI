#!/usr/bin/perl
#Helps format the fasta (nucleotide) file in a way that is compatible with the remote blast search script
#Expects the file to be named strain.fasta and the input of a genome with -genome (Genome) is required
#OVERWRITES OLD FILE. MAKE SURE YOU KEEP A BACKUP!!!!!!
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Spec;

my $help   = 0;
my $man    = 0;
my $genome = '';
my $o      = '';
my $single = '';
my $split  = 0;
my $scri   = 0;

my $signal = GetOptions(
                         'help'        => \$help,
                         'man'         => \$man,
                         'genome=s'    => \$genome,
                         'overwrite|o' => \$o,
                         'single|S'    => \$single,
                         'split'       => \$split,
                       );

my @infiles = @ARGV;
my $count   = @infiles;

#Print help statements if specified
pod2usage( -verbose => 1,
           -output  => *STDERR )
  if $help == 1;
pod2usage( -verbose => 2,
           -output  => *STDERR )
  if $man == 1;

pod2usage(
     -verbose => 1,
     -output  => *STDERR,
     -msg =>
       "No genome specified. Please specify genome with -genome (genomename).\n"
  )
  if !$genome && !$split && !$scri;
pod2usage(
    -verbose => 1,
    -output  => *STDERR,
    -msg =>
      "No input files specified. Please specify files using $0 -genome (genomename) infile1.fasta infile2.fasta etc.\n"
  )
  if ( scalar(@ARGV) == 0 );

pod2usage(
    -verbose => 0,
    -output  => *STDERR,
    -msg =>
      "You must provide the -overwrite command to run the script.  This ensures you are ready to overwrite your old files.  Make sure you make a backup of your old file before generating a new fasta file with this script.  THIS IS ESSENTIAL OR YOU WILL LOSE YOUR DATA!!!\n\nAlternatively, use the -single flag to only write one file. These options are exclusive to one another, with -overwrite taking precedence."
  )
  if ( ( !$o && !$single || $o && $single ) && !$split );

if ($o) {
    print STDERR "Overwrite flag given.\n";
    print STDERR "Attempting to rename $count files\n";
} elsif ($single) {
    die("Multiple input files given, with the single flag specified.  Retry with only one input file.\n"
       )
      if $count != 1;
    print STDERR "Single flag given.  Printing output to STDOUT\n";
}

my $strain = '';

if ( $o || $single ) {
    &rename;
}

if ($split) {
    &split;
}

sub rename {

    foreach my $infile (@infiles) {
        open INFILE, "$infile" or die "$infile unavailable : $!";

        my $fh;
        if ($o) {
            open( $fh, '>', "temp" ) or die "Unable to open temp : $!";
        } else {
            $fh = \*STDOUT;
        }

        my ($vol, $dir, $file) = File::Spec->splitpath( $infile );

        $strain = $file;
        $strain =~ s/\.[^\.]+$//;
        my $temp = $strain;
        $temp =~ s/_/ /g;
        while (<INFILE>) {
            my $line = $_;
            chomp($line);
            if ( $line =~ /^>/ ) {
                my @data = split( '\|', $line );
                if ( scalar(@data) > 1 ) {
                    for ( my $i = 0 ; $i < scalar(@data) ; $i++ ) {
                        $data[$i] =~ s/>//;
                        if ( $i == 0 ) {
                            print $fh '>gnl|' . $strain . '|';
                        } elsif ( $i == ( scalar(@data) - 1 ) ) {
                            print $fh $data[$i] . " [$genome $temp]\n";
                            last;
                        }

                        print $fh $data[$i] . '|';
                    }
                } else {
                    $line = substr $line, 1;
                    print $fh ">gnl|$strain|$line [$genome $temp]\n";

                }
            } else {
                print $fh $_;
            }
        }
        close $fh;
        system("rm -f $infile")   if $o;
        system("mv temp $infile") if $o;
    }

}

sub split {

    foreach my $infile (@infiles) {
        if ( !-e $infile ) {
            print STDERR "Cannot find $infile.\n Skipping...\n";
            next;
        }
        open INFILE, "$infile" or die "$infile unavailable : $!";
        my %output;
        my $i = 0;
        while (<INFILE>) {
            my $line = $_;
            chomp($line);
            if ( $line =~ /^>/ ) {
                $i++;
                my @data = split( '\|', $line );
                if ( $data[1] ) {
                    $data[1] .= "_$i";
                    $strain = $data[1];
                } else {
                    $strain = $line . "_$i";
                }
                $output{$strain} = [];
                $line = join( '|', @data );
                push( @{ $output{$strain} }, $line );
            } else {
                push( @{ $output{$strain} }, $line );
            }
        }
        close INFILE;
        foreach my $strain ( keys %output ) {
            open OUTPUT, ">$strain.fna"
              or die "$strain.fna is unavailable : $!";
            print OUTPUT join( "\n", @{ $output{$strain} } );
            close OUTPUT;
        }
    }

}

__END__

=head1

=head1 NAME

fasta_format.pl - Renames nucleotide fasta files with proper names for blast searching and MLSA - B<OVERWRITES OLD FILE!>

=head1 SYNOPSIS

fasta_format.pl (-overwrite or -single) -genome 'Genus species' strain1.fasta strain2.fasta

Expects the strain name to be included in the filename (e.g. K12.fasta) and the genus, species name provided using -genome (e.g. 'Escherichia coli' - must be in quotes)

=head1 OPTIONS

=over 8

=item B<-help|h>

Print a brief help message and exits.

=item B<-man>

Print a verbose help message and exits.

=item B<-genome> 'Genus species'

Provide an organism name to append to the strain name provided by the input file.

=item B<-overwrite> or B<-single>

-overwrite is required for generating multiple new fasta files at a time.  -single can be used to perform one rename at at time, with output to STDOUT instead of a file.

=back

=head1 DESCRIPTION

This program will rename the header line in a provided fasta file such that the data can be used in the MLSA pipeline.  The MLSA pipeline assumes, as with NCBI, the format of 

>gi|12345679|otherID| description [Genus species strain]

This program will generate a file following that format. eg :

>gnl|K12| chromosome 1 [Escherichia coli K12]

The only important parts for the MLSA script are those found at the start of the header.  The script looks for lines starting with >gnl| and takes the strain name after the first pipe and then takes the species name from inside the brackets.  Can handle subspecies as well (Salmonella enterica subsp. enterica) if those are given as the -genome name.

=cut
