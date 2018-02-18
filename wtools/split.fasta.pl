#!/usr/bin/env perl

# split.fasta.pl
# 
#
# Author : mieg
#
# Date : octobre 1999
#

# Takes a big fasta file and exports in dir ./Dna
#  one fasta file per sequence
#ENVIRONMENT

$| = 1; #turns off output buffering for STDOUT.

use Getopt::Std;

# get command line arguments
&getopts('f:');
unless ($opt_f) {
    # must have sequence argument
    print "Usage:  split.fasta.pl -f fasta_file\n".
    "Takes a big fasta file and exports in dir ./Dna\n".
    "one fasta file per sequence\n" ;
    exit 1;
}

open (DNA, "/bin/ls Dna |") || die "cannot find subdirectory ./Dna" ;
close DNA ;

open (FASTA, $opt_f)  || die "cannot open $opt_f" ;
while (<FASTA>)
{
    chop ;

    if (/^\>(.*)/)
	{
	    $name = $1 ;
	    open (SEQ, "> Dna/$name") || die "cannot open Dna/$name" ; 
	    print SEQ ">$name\n" ;
	}
    else
    {
	print SEQ "$_\n";
    }

}
close SEQ ;



