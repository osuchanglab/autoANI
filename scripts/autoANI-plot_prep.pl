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

my $headers = 0;

GetOptions("headers" => \$headers);

my $infile = shift;
my $outfile = $infile.".divergence";
my $names = $infile.".names";

open INFILE, "$infile" or die "$infile unavailable : $!";

open my $outfh, ">", "$outfile" or die "Unable to open $outfile for writing: $!";

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
            $data = 0;
        } else {
            $data = 100 - $data;
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
