#!/usr/bin/perl
use strict;
use warnings;
use Term::ANSIColor qw(:constants);
   $Term::ANSIColor::AUTORESET=1;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/";
use FASTA;
use ANNOTATION;

my $usage=<<USAGE;
	Usage:	$0 <input.gff> [-o <output.gff>]
USAGE

print GREEN join ' ', ($0,@ARGV,"\n");

############################################################

my ($oput_gff);
GetOptions ("o:s"=>\$oput_gff);
die $usage unless @ARGV;
my (@iput) = @ARGV;
$oput_gff ||= 'out.gff';


my %gff;
foreach (@iput) {
	ANNOTATION::Read_GFF ($_, \%gff);
}


my $count;
open my $OPUT, ">$oput_gff" or die;
foreach my $scaffold_id (ANNOTATION::sort_scaffold_id(\%gff)) {
	my $gff = $gff{$scaffold_id};
	best_gene ($gff);
	ANNOTATION::Print_GFF ($OPUT, $gff, $scaffold_id);
}
close $OPUT;


sub best_gene {
	my $gff = shift;
	my %len;
	foreach my $g (@$gff) {
		my $len = 0;
		my $gene_id = $$g[0]{id};
		my $mRNA_id = $$g[1]{id};
		foreach (@$g) {
			#print Dumper($_);
			$len += $$_{end} - $$_{start} + 1 if ($$_{type} eq 'CDS');
		}
		$len{$gene_id}{len} ||= 0;
		if ($len > $len{$gene_id}{len}) {
			$len{$gene_id}{len} = $len;
			$len{$gene_id}{mRNA_id} = $mRNA_id; 
		}
	}
	foreach my $g (@$gff) {
		my $gene_id = $$g[0]{id};
		my $mRNA_id = $$g[1]{id};
		if ($mRNA_id ne $len{$gene_id}{mRNA_id}) {
			$$g[0]{filter} = 1;
		}
		else {
			$$g[0]{type} = 'gene';
			$len{$gene_id}{mRNA_id} = '';
		}
	}
}

__END__
