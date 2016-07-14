#!/usr/bin/perl
#######################################################################
#
# COPYRIGHT NOTICE
#
# autoANI-blast.pl - script to automate BLAST calls from autoANI.pl
# Splitting the script in this way makes parallelization a bit easier
# from my perspective
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

my $coverage   = shift;
my $size       = shift;
my $pid_cutoff = shift;
my $output     = shift;
my @command    = join( " ", @ARGV );
my %results;
my %genomes;

@command = join( " ", @command, '|', 'tee', $output );

my @output = `@command`;

#foreach my $line (@output) {
#    chomp($line);
#    next if ( $line =~ /#/ );
#    my @data      = split( "\t", $line );
#    my $query     = $data[0];
#    my $subject   = $data[1];
#    my $percentid = $data[2];
#    my $length    = $data[3];
#    $query =~ s/_[0-9]+(?=$)//;
#
#    my $qacc = get_accn($query);
#    my $sacc = get_accn($subject);
#
#    $qacc =~ s/\.[0-9]+//;
#    $sacc =~ s/\.[0-9]+//;
#
#    if ( $query ne $subject ) {
#        next
#          if (    $length < ( ( $coverage / 100 ) * $size )
#               || $percentid < $pid_cutoff );
#        push( @{ $results{$qacc}{$sacc} }, $percentid );
#    }
#
#}
#
#
#open my $tmpout, ">>", "ani.tmp" or die "Unable to open ani.tmp for writing!\n";
#
#foreach my $genome ( sort keys %results ) {
#    foreach my $subject ( sort keys %{ $results{$genome} } ) {
#        my ( $ani, $hits );
#        my $sum  = 0;
#        my @data = @{ $results{$genome}{$subject} };
#        if ( !@data ) {
#            $ani  = 'N/A';
#            $hits = 0;
#        } else {
#            $hits = @data;
#            $sum += $_ for @data;
#            $ani = $sum / $hits;
#        }
#        print $tmpout
#          join( "\t", $genome, $subject, sprintf( "%.3f", $ani ), $hits )
#          . "\n";
#    }
#}
#
#close $tmpout;
#
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
