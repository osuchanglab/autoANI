#!/usr/bin/env perl

use warnings;
use strict;
use v5.10;
use File::Spec;
use Net::FTP;
use Env qw(PATH);

if ( -d 'edirect') {
    say "edirect folder already found!";
    say "Downloading and updating edirect files.";
} else {
    say "Downloading and installing edirect files.";
}

if ( -e 'edirect.zip' ) {
    system('rm -f edirect.zip');
}

my $ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1);
$ftp->login;
$ftp->binary;
$ftp->get("/entrez/entrezdirect/edirect.zip");

if ( -e 'edirect.zip' ) {
    say "Successfully downloaded edirect.zip";
    say "Unzipping edirect.zip";
} else {
    say "Unable to download edirect.zip.";
    say "Check your permissions and try again.";
    exit(-1);
}

system('unzip -u edirect.zip') == 0 or die "Unable to unzip edirect.zip : $!";

system('rm -f edirect.zip');

`sed -i 's/perl -w/env perl\\nuse warnings;/' edirect/ftp-cp`;
`sed -i 's|/usr/bin/perl|/usr/bin/env perl|' edirect/edirect.pl`;

say "Running ./edirect/setup.sh file!";
say "--------------------------------";

system('./edirect/setup.sh');

say "\n--------------------------------";
say "Finished running setup.sh file!\n";

my $edir = File::Spec->rel2abs('./edirect');

my $edir_search = qr!$edir!;

if ( $PATH !~ $edir_search ) {
    say "Add $edir to your path environment variable to use these E-Utilities (see above for command for Bash).";
    say "On tcsh or csh, add this to your ~/.cshrc file:";
    say "setenv PATH \${PATH}:$edir";
} else {
    say "$edir found in PATH. You are ready to use edirect E-Utilities!";
}

