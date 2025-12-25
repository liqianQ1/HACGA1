#!/usr/bin/perl
use strict;
use warnings;
use Term::ANSIColor qw(:constants);
   $Term::ANSIColor::AUTORESET=1;
use Getopt::Long;

use FindBin qw($Bin);
use lib "$Bin/";
use FASTA;
use ANNOTATION;

my $usage=<<USAGE;
	Usage:	$0 <pasa.orfs.gff> [-o <output.gff>] [-p <input.pep>] [-c <input.cds>] [-b input.blast] [-min <min_length>] [-max <max_length>] [-cds <int>]
USAGE

print GREEN join ' ', ($0,@ARGV,"\n");

############################################################

my ($oput_gff, $pep, $cds, $blast, $min_len, $max_len, $cds_num);
GetOptions ("o:s"=>\$oput_gff,
            "p:s"=>\$pep,
            "c:s"=>\$cds,
            "b:s"=>\$blast,
            "min:i"=>\$min_len,
            "max:i"=>\$max_len,
            "cds:i"=>\$cds_num);
die $usage unless (@ARGV);
my ($pasa) = @ARGV;
$oput_gff ||= 'gene.filter.gff';
$min_len  ||= 100;
$max_len  ||= 99999;
$cds_num  ||= 1;


my $oput = $oput_gff;
$oput =~ s/\.\w+(.tmp)?$//;
my $oput_cds = "$oput.cds";
my $oput_pep = "$oput.pep";


my %gff;
ANNOTATION::Read_GFF ($pasa, \%gff);


my %blast;
if ($blast) {
	my @blast;
	ANNOTATION::Read_Blast ($blast, \@blast);
	foreach (@blast) {
		my $qlen = $$_{Query_end} - $$_{Query_start} + 1;
		my $slen = $$_{Subject_end} - $$_{Subject_start} + 1;
		if ($qlen==$$_{Query_len} && $slen==$$_{Subject_len} && $$_{Query_len}>=$min_len && $$_{Query_len}<=$max_len) {
			my $id = $$_{Query_ID};
			$blast{$id} = 1;
		}
	}
}


my %pep;
if ($pep) {
	open my $IPUT, "<$pep" or die;
	while (!eof $IPUT) {
		my $fa = FASTA::read_fasta_one_seq ($IPUT);
		$pep{$$fa{id}} = $fa;
	}
	close $IPUT;
}


my %cds;
if ($cds) {
	open my $IPUT, "<$cds" or die;
	while (!eof $IPUT) {
		my $fa = FASTA::read_fasta_one_seq ($IPUT);
		$cds{$$fa{id}} = $fa;
	}
	close $IPUT;
}


foreach my $gff (values %gff) {
	foreach (@$gff) {
		my $cds = 0;
		my $cds_len = 0;
		my $utr5 = 0;
		my $utr3 = 0;
		foreach (@$_) {
			if ($$_{type} eq 'CDS') {
				$cds++;
				$cds_len += $$_{end} - $$_{start} + 1;
			}
			elsif ($$_{type} eq 'five_prime_UTR') {
				$utr5++;
			}
			elsif ($$_{type} eq 'three_prime_UTR') {
				$utr3++;
			}
		}
		$$_[0]{filter} = 1 if ($cds<$cds_num or ! $utr3 or !$utr5 or $cds_len<$min_len*3+3 or $cds_len>$max_len*3+3);
	}
}


my %count;
open my $OGFF, ">$oput_gff" or die;
open my $OCDS, ">$oput_cds" or die if ($cds);
open my $OPEP, ">$oput_pep" or die if ($pep);
foreach my $scaffold_id (ANNOTATION::sort_scaffold_id(\%gff)) {
	my $gff = $gff{$scaffold_id};
	ANNOTATION::sort_gff ($gff);
	ANNOTATION::sort_gene ($gff);
	foreach (@$gff) {
		my $id = $$_[1]{id};
		$$_[0]{filter} = 1 if ($blast and !$blast{$id});
		if (!$$_[0]{filter}) {
			FASTA::write_fasta_one_seq ($OCDS, $cds{$id}) if ($cds);
			FASTA::write_fasta_one_seq ($OPEP, $pep{$id}) if ($pep);
		}
	}
	ANNOTATION::Print_GFF ($OGFF, $gff, $scaffold_id);
}
close $OGFF;
close $OCDS if ($cds);
close $OPEP if ($pep);

__END__
