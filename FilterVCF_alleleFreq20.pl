#!/usr/bin/perl -w
#
#	This script filters SNP's in a VCF file to only retain
#	SNP's with a frequency between 0.2 and 0.8. It uses the
#	AF (allele frequency of alternate allele) info for filtering
#
#	Usage:
#	perl FilterVCF_alleleFreq20.pl VCF_filename > Outputfile
#
#	This script was written by Henrik H. De Fine Licht
#
################################################################
use strict;

my $freq;
my $filename;

$filename = $ARGV[0];

open(my $fh, '>', $filename) or die "Couldn't open file $filename, $!";

while (my $line = <>) {
    if ( $line =~ /AF1=(\d\.?\d*);/ ) {
	$freq = $1;
	if ($freq < 0.8 and $freq > 0.2 ) {
	    print $line;
	    print $fh $freq,"\n";
		}
    }
}

close $fh;
