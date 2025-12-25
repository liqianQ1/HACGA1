#!/usr/bin/env perl
use strict;
use warnings;
use Term::ANSIColor qw(:constants);
   $Term::ANSIColor::AUTORESET=1;
use Getopt::Long;

#use FindBin qw($Bin);
#use lib "$Bin/../lib/";
#use FASTA;
#use ANNOTATION;

my $usage=<<USAGE;
	Usage:	$0 <augustus.gff> [-o <output.gff>]
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
	Read_GFF ($_, \%gff);
}


my $count;
open my $OPUT, ">$oput_gff" or die;
foreach my $scaffold_id (sort_scaffold_id(\%gff)) {
	my $gff = $gff{$scaffold_id};
	sort_gff ($gff);
	motify ($_) foreach (@$gff);
	unique_gene_id ($gff, \$count, $scaffold_id);
	Print_GFF ($OPUT, $gff, $scaffold_id);
}
close $OPUT;


sub motify {
	my $g = shift;
	my @g;
	foreach (my $i=0;$i<@$g;$i++) {
		if ($$g[$i]{type} eq 'CDS' && $$g[$i-1]{type} ne 'exon') {
			my $j = @g;
			$g[$j]{program} = $$g[$i]{program};
			$g[$j]{type}    = 'exon';
			$g[$j]{start}   = $$g[$i]{start};
			$g[$j]{end}     = $$g[$i]{end};
			$g[$j]{score}   = '.';
			$g[$j]{strand}  = $$g[$i]{strand};
			$g[$j]{phase}   = '.';
			$g[$j]{parent}  = $$g[$i]{parent};
		}
		push @g, $$g[$i] if ($$g[$i]{type} ne 'intron');
	}
	$$_{ann} = '' foreach (@g);
	@$g = @g;
}


############################################################

sub sort_scaffold_id {
	my $scaffold = shift;
	my %scaffold_id;
	my $flag = 100000;
	foreach (keys %$scaffold) {
		#/(\d+)/ or die;
		if(/(\d+)/){
			$scaffold_id{$_} = $1;
		}else{
			$scaffold_id{$_} = $flag;
			$flag++;
		}
	}
	my @scaffold_id = sort {$scaffold_id{$a}<=>$scaffold_id{$b}} keys %scaffold_id;
	return @scaffold_id;
}


sub sort_gff {
	my $gff = shift;
	@$gff = sort {$$a[0]{start}<=>$$b[0]{start} || $$b[0]{end}<=>$$a[0]{end} || $$a[0]{score}=~/\d/ && $$b[0]{score}=~/\d/ && $$b[0]{score}<=>$$a[0]{score}} @$gff;
}

############################################################

sub Read_GFF {
	my ($iput, $gff) = @_;
	open IPUT, "<$iput" or die;
	while (<IPUT>) {
		/^\w/ or next;
		chomp;
		my @c = split /\t/;
		my $scaffold_id = $c[0];
		my %h;
		$h{program} = $c[1];
		$h{type}    = $c[2];
		$h{start}   = $c[3];
		$h{end}     = $c[4];
		$h{score}   = $c[5];
		$h{strand}  = $c[6];
		$h{phase}   = $c[7];
		$h{ann}     = $c[8]||'';
		($h{start}, $h{end}) = ($h{end}, $h{start}) if ($h{start}>$h{end});
		if ($h{ann} =~ /^ID=([^;]+);?(.*)$/) {
			$h{id}  = $1;
			$h{ann} = $2||'';
		}
		if ($h{ann} =~ /^Parent=([^;]+);?(.*)$/) {
			$h{parent} = $1;
			$h{ann}    = $2||'';
		}

		if (!$h{parent}) {
			my $i = $#{$$gff{$scaffold_id}} + 1;
			$$gff{$scaffold_id}[$i][0] = \%h;
		}
		elsif ($h{parent} eq $$gff{$scaffold_id}[-1][0]{id} && $#{$$gff{$scaffold_id}[-1]}>1 && $h{type} eq 'mRNA') {
			my $i = $#{$$gff{$scaffold_id}} + 1;
			%{$$gff{$scaffold_id}[$i][0]} = %{$$gff{$scaffold_id}[$i-1][0]};
			$$gff{$scaffold_id}[$i][0]{type} = '';
			$$gff{$scaffold_id}[$i][1] = \%h;
		}
		else {
			push @{$$gff{$scaffold_id}[-1]}, \%h;
		}
	}
	close IPUT;
}


sub unique_gene_id {
	my ($gff, $count, $scaffold_id, $program) = @_;
	$program ||= $$gff[0][0]{program};
	my $number;
	my $parent;
	my $m = 1;
	foreach my $g (@$gff) {
		$$g[0]{filter} and next;
		if ($$g[0]{type} eq 'gene') {
			$number = sprintf "%05d", ++$$count;
			$$g[0]{program} = $program;
			$$g[0]{id} = "$scaffold_id.g$number";
			$parent = $$g[0]{id};
			$m = 1;
		}
		elsif (!$$g[0]{type}) {
			$$g[0]{program} = $program;
			$$g[0]{id} = "$scaffold_id.g$number";
			$parent = $$g[0]{id};
		}
		else {
			die;
		}
		my %n;
		$n{exon} = 0;
		$n{utr5p} = 0;
		$n{utr3p} = 0;
		for (my $i=1;$i<@$g;$i++) {
			my $type = $$g[$i]{type};
			if ($type eq 'mRNA') {
				$$g[$i]{program} = $program;
				$$g[$i]{id} = "$parent.m$m";
				$$g[$i]{parent}  = $$g[0]{id};
				$parent = $$g[$i]{id};
				$m++;
				$n{$_} = 0 foreach (keys %n);
			}
			else {
				if ($type eq 'CDS') {
					$type = 'cds';
				}
				elsif ($type eq 'five_prime_UTR') {
					$type = 'utr5p';
				}
				elsif ($type eq 'three_prime_UTR') {
					$type = 'utr3p';
				}
				$$g[$i]{program} = $program;
				$$g[$i]{parent} = $parent;
				$$g[$i]{id} = "$parent.$type";
				if ($type =~ /^exon$|^utr5p$|^utr3p$/) {
					$n{$type}++;
					$$g[$i]{id} .= $n{$type};
				}
			}
		}
	}
}

sub Print_GFF {
	my ($OPUT, $gff, $scaffold_id, $scaffold_len) = @_;
	print $OPUT "##sequence-region $scaffold_id 1 $scaffold_len\n" if ($scaffold_len);
	foreach my $g (@$gff) {
		$$g[0]{filter} and next;
		foreach (@$g) {
			$$_{type} or next;
			my $txt = '';
			$txt .= "$scaffold_id\t$$_{program}\t$$_{type}\t$$_{start}\t$$_{end}\t$$_{score}\t$$_{strand}\t$$_{phase}\t";
			$txt .= "ID=$$_{id};" if (defined $$_{id});
			$txt .= "Parent=$$_{parent};" if (defined $$_{parent});
			$txt .= "$$_{ann}" if ($$_{ann});
			print $OPUT "$txt\n";
		}
	}
}

__END__
