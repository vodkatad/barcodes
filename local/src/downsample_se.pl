#!/usr/bin/perl
use strict;
use warnings;
#use PerlIO::gzip;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
##
## Usage ./sampleFastq.pl <fastq r1> <fastq r2> <outFastq r1> <outFastq r2> <prob of keeping reads>
## modified from https://www.biostars.org/p/110107/
#
if (scalar(@ARGV) != 3) {
	die "Usage: ./sampleFastq.pl <fastq r1> <outFastq r1> <prob of keeping reads>";
}

srand(42);
my $fqf = new IO::Uncompress::Gunzip $ARGV[0], MultiStream => 1 or die "could not read $ARGV[0]: $!";
my $fqoutf = new IO::Compress::Gzip $ARGV[1] or die "could not write $ARGV[2]: $!";
my $proba = $ARGV[2];

print STDERR "Keeping $proba reads\n";

my $nbLines = 1;
my $fqRecordf = '';
while (my $line1=<$fqf> ){
    #print STDERR "$nbLines $line1\n";
    $fqRecordf .= $line1;
    if ($nbLines % 4 == 0) {
        my $random = rand(1);
        if ($random <= $proba) {
            print $fqoutf $fqRecordf;
        }
        $fqRecordf = '';
    }
    $nbLines++;
}
close $fqf;
close $fqoutf;
