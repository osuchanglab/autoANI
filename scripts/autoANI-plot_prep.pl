#!/usr/bin/env perl
#######################################################################
#
# COPYRIGHT NOTICE
#
# autoANI-plot_prep.pl - A script to automate ANI divergence (100-ANI)
# calculations.  Generates a matrix file (.divergence) and names file
# (.names) suitable for input to R or other software for heat map 
# generation and hierarchical clustering.
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

my $headers = 0;
my $divergence = 0;
my $help = 0;

GetOptions(
    "headers" => \$headers,
    "d|divergence" => \$divergence,
    "h|help" => \$help
);

if ($help) {
    pod2usage(
        -verbose => 1,
        -output => \*STDERR
    );
}

my $infile = shift;
#my $outfile = $infile.".divergence";
my $names = $infile.".names";

open INFILE, "$infile" or die "$infile unavailable : $!";

#open my $outfh, ">", "$outfile" or die "Unable to open $outfile for writing: $!";
my $outfh = \*STDOUT;

open my $namesfh, ">", "$names" or die "Unable to open $names for writing: $!";

my $headerline = <INFILE>;

if ($headers == 1) {
    print $outfh $headerline;
}

my @names;

while (<INFILE>) {
    my $line = $_;
    chomp($line);
    my ($name, @data) = split( "\t", $line );
    print $namesfh $name."\n";
    foreach my $data (@data) {
        if ( $data eq '----------' ) {
            $data = 100;
            if ( $divergence ) {
                $data = 0;
            }
        } else {
            if ( $divergence ) {
                $data = 100 - $data;
            }
        }
    }
    if ( $headers == 1 ) {
        print $outfh $name."\t";
    }
    print $outfh join( "\t", @data ) . "\n";
}
close INFILE;

close $outfh;

close $namesfh;

__END__

=head1 Name

autoANI-plot_prep.pl - Help printing matrices for plotting/heatmap generation in R

=head1 SYNOPSIS

autoANI-plot_prep.pl [options] ani.out > ani_plot.out

=head1 OPTIONS

=item B<-help>

Print a brief help message and exits.

=item B<-divergence>

Print ANI divergence (100-ANI) instead of ANI.

=item B<-headers>

Print header lines. By default, no headers are printed.

