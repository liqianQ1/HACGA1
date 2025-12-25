#!/usr/bin/perl -w
use warnings;
use strict;
use Bio::SeqIO;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

my $in = Bio::SeqIO->new(-file    =>    "$infile",
			 -format  =>    "fastq");
my $out = Bio::SeqIO->new(-file    =>   ">$outfile",
			  -format  =>    "fasta");

while (my $seq = $in->next_seq()){
	$out->write_seq($seq);
}

