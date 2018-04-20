#!/usr/bin/env perl
#######################################################################
#
# COPYRIGHT NOTICE
#
# autoANI-sort.pl - A script to sort the output from autoANI.pl
# This allows the user to control the reference genome in the output
# and sort all of the values in a 2 dimensional matrix.  The generated
# output can then be visualized using spreadsheet software in a
# meaningful way.
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
use Getopt::Long;
use Pod::Usage;

my $infile;
my $hitsfile;
my $ref = '';
my $help = 0;
my %notseen;
my %hits;
my %values;
my $hitlimit = 0;
my @hitfilter;
my $sample;
my %keep;

GetOptions( 'reference=s' => \$ref,
            'help|h'        => \$help,
            'hitsfile=s'      => \$hitsfile,
            'hitlimit=i'      => \$hitlimit,
            'subsample=s' => \$sample);

if ($help) {
    pod2usage(
        -verbose => 1,
        -output => \*STDERR );
}

$infile = shift;

if ( !$ref && !$infile ) {
    print STDERR
      "No reference genome provided! Please provide a reference genome to sort by with -reference [Genome Name]\n";
      exit(-1);
}

if ( !-e $infile ) {
    print STDERR
      "Unable to find input file $infile. Check your input and try again.\n";
    exit(-1);
}

if ($hitsfile) {
    if ( !-e $hitsfile ) {
        print STDERR
          "Unable to file hits file $hitsfile.  Check input and try again.\n";
        exit(-1);
    }

    open my $hitsfh, "<", $hitsfile
      or die "Unable to open hits file $hitsfile : $!";

    getData( $hitsfh, \%hits );

    close $hitsfh;

}

if ($sample) {
    if ( ! -e $sample ) {
        print STDERR "Please provide a valid, existing file name for -sample argument and try again\n";
        exit(-1);
    }

    open my $samplefh, "<", $sample or die "Unable to open sample file $sample : $!";
    while ( <$samplefh> ) {
        my $line = $_;
        chomp($line);
        $keep{$line} = 1;
    }
}

open my $infh, "<", $infile or die "Unable to open input file $infile : $!";

getData($infh,\%values);

close $infh;

if ( ! $ref || ! exists $values{$ref} ) {
    my @refs = sort keys %values;
    print STDERR "Unable to find reference name \"".$ref."\" in data.\n";
    print STDERR "Possible values are:\n".join("\n",@refs)."\n";
    exit(-1);
}

foreach my $genome ( keys %hits ) {
    my $n = 0;
    my $total = (keys %hits) - 1;
    foreach my $subject ( keys %hits) {
        next if $genome eq $subject;
        $n += $hits{$genome}{$subject};
    }
    my $avg = $n/$total;
    if ( $avg <= $hitlimit) {
        push(@hitfilter,$genome);
    }
}

if (scalar(@hitfilter) > 0) {

    print STDERR "These genomes are filtered due to low number of average hits:\n";
    print STDERR join("\n",@hitfilter)."\n";

    foreach my $filtered (@hitfilter) {
        delete $values{$ref}{$filtered};
    }

}

my @neworder;
push(@neworder, $ref);
#my @sortorder =  sort { $values{$ref}{$b} <=> $values{$ref}{$a} } keys %{$values{$ref}};
#my $sortref = shift @sortorder;
%notseen = map { $_ => 1 } sort keys %{$values{$ref}};
#@sortorder = ();
while ( %notseen ) {
   my @sortorder = sort { $values{$ref}{$b} <=> $values{$ref}{$a} } keys %notseen;
   $ref = shift (@sortorder);
   push(@neworder, $ref);
   delete($notseen{$ref});
}

print join("\t","",@neworder)."\n";
foreach my $genome (@neworder) {
    print $genome;
    foreach my $genome2 (@neworder) {
        if ($genome eq $genome2) {
            print "\t"."100";
        } else {
            if ( exists($values{$genome}{$genome2}) ) {
                print "\t".$values{$genome}{$genome2};
            } else {
                print STDERR "Unable to find data for $genome2 in $genome\n";
                print STDERR "Check your input file for missing values and try again\n";
                exit;
            }
        }
    }
    print "\n";
}

sub getData {
    my $fh   = shift;
    my $hash = shift;
    my $i   = 0;
    my $row = 0;
    my @order;

    while (<$fh>) {
        my $line = $_;
        chomp($line);
        if ( $i == 0 ) {
            my $junk;
            ( $junk, @order ) = split( "\t", $line );
        } else {
            my $line = $_;
            chomp($line);
            my ( $name, @data ) = split( "\t", $line );
            if ( scalar(@data) != scalar(@order) ) {
                print STDERR "Missing data point for $name. Check to make sure all of the data is in this row and is formatted correctly\n";
                exit;
            }
            for ( my $j = 0 ; $j < scalar(@data) ; $j++ ) {
                if ( $name eq $order[$j] ) {
                    next;
                } elsif ( %keep ) {
                    next if ! exists( $keep{$name} );
                    next if ! exists( $keep{$order[$j]} );
                }
                $hash->{$name}{ $order[$j] } = $data[$j];
            }
        }
        $i++;
    }
}

__END__

=head1 Name

autoANI-sort.pl - Sort ANI results table by supplied reference genome

=head1 SYNOPSIS

autoANI-sort.pl -reference "Genus species" input.tab

=head1 OPTIONS

If no reference is supplied, a list of possible values is printed.

=item B<-help>

Print a brief help message and exits.

=item B<-reference>

Genome name, in quotes, to sort by.


