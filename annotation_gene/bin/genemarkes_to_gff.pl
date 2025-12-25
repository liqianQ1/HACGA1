#!/usr/bin/perl -w
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
	Usage:	$0 <scaffold.masked.fa> <genemarkes.gff> [-o <output.gff>]
USAGE

print GREEN join ' ', ($0,@ARGV,"\n");

############################################################

my ($oput_gff);
GetOptions ("o:s"=>\$oput_gff);
die $usage unless (@ARGV >= 2);
my ($iput, @iput) = @ARGV;
$oput_gff ||= 'out.gff';


my %gff;
foreach (@iput) {
	Read_GFF ($_, \%gff);
}


my $count;
open my $OPUT, ">$oput_gff" or die;
open my $IPUT, "<$iput" or die;
while (!eof $IPUT) {
	my $scaffold = read_fasta_one_seq ($IPUT);
	my $gff = $gff{$$scaffold{id}};
	$gff or next;
	motify ($_) foreach (@$gff);
	sort_gff ($gff);
	sort_gene ($gff);
	extract_fa ($gff, $scaffold);
#	filter_gap ($gff, 1000);
	unique_gene_id ($gff, \$count, $$scaffold{id}, 'GeneMark.hmm');
	Print_GFF ($OPUT, $gff, $$scaffold{id});
}
close $OPUT;
close $IPUT;


sub motify {
	my $g = shift;
	$$_{ann} =~ s/Name=[^;]+\;?// foreach (@$g);
}

sub read_fasta_one_seq {
		my $IPUT = shift;
		my $fasta;
		$/="\n";
		<$IPUT> =~ /\>?(.*?)\s+(.*)?$/ or die;
		$$fasta{id}  = $1;
		$$fasta{ann} = $2||'';
		$/='>';
		$$fasta{seq} = <$IPUT>;
		chomp $$fasta{seq};
		$$fasta{seq} =~ s/\n//g;
		$$fasta{seq} =~ s/\*$//;
		$$fasta{len} = length $$fasta{seq};
		$/="\n";
		return $fasta;
}


sub sort_gff {
	my $gff = shift;
	@$gff = sort {$$a[0]{start}<=>$$b[0]{start} || $$b[0]{end}<=>$$a[0]{end} || $$a[0]{score}=~/\d/ && $$b[0]{score}=~/\d/ && $$b[0]{score}<=>$$a[0]{score}} @$gff;
}


sub sort_gene {
	my $gff = shift;
	foreach my $g (@$gff) {
		@$g = sort {$$a{start}<=>$$b{start} || $$b{end}<=>$$a{end}} @$g;
	}
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




sub Read_rRNA2{
	my ($iput, $gff) = @_;
	open FH,"cat $iput | tr -s \" \" | sed  \"s\/ \/\t\/g\" | " 
	|| die "Can't open file";
	my $cmfetch = "/public/software/Infernal/bin/cmfetch";
	my $Rfam = "/public/database/Rfam/CMs/Rfam.cm";
	while (<FH>) {
		next unless !/^#/;
		chomp;
		my %h;
		my $coverage;
		my ($target_name,$Taccession,$query_name,
		$Qaccession,$mdl,$mdl_from,$mdl_to,
		$seq_from,$seq_to,$strand,$trunc,
		$pass,$gc,$bias,$score,$Evalue,$inc,$description_of_target)  = split "\t";
		my $target_len = `$cmfetch $Rfam  $Taccession | grep CLEN | tr -s ' ' | cut -d ' ' -f2`;
		($h{strand}, $h{start}, $h{end}, $h{score}) = ($strand , $seq_from , $seq_to ,$score);

		$coverage = sprintf "%.2f", 100*($mdl_to-$mdl_from+1)/$target_len; ###
		$h{program} = 'INFERNAL';
		$h{type}    = 'rRNA';
		$h{phase}   = '.';
		$h{ann}     = "Target=$target_name $mdl_from $mdl_to;coverage=$coverage;";
		my $i = $#{$$gff{$query_name}} + 1;
		$$gff{$query_name}[$i][0] = \%h;
	}
	close FH;
}

sub Read_rRNA {
	my ($iput, $gff) = @_;
	open IPUT, "<$iput" or die;
	while (<IPUT>) {
		/^#/ and next;
		chomp;
		my $scaffold_id;
		my %h;
		my ($target, $target_len, $target_start, $target_end, $coverage);
		($target, $target_len, $target_start, $target_end, $h{strand}, $scaffold_id, $h{start}, $h{end}, $h{score}) = (split /\t/, $_)[0,1,2,3,4,5,7,8,10];
		$target =~ /(rRNA_[\d\.]+S)/ or die;
		$target = $1;
		$coverage = sprintf "%.2f", 100*($target_end-$target_start+1)/$target_len;
		$h{program} = 'BLASTN';
		$h{type}    = 'rRNA';
		$h{phase}   = '.';
		$h{ann}     = "Target=$target $target_start $target_end;coverage=$coverage;";
		my $i = $#{$$gff{$scaffold_id}} + 1;
		$$gff{$scaffold_id}[$i][0] = \%h;
	}
	close IPUT;
}


sub Read_Solar {
	my ($iput, $gff) = @_;
	open IPUT, "<$iput" or die;
	while (<IPUT>) {
		/^#/ and next;
		chomp;
		my ($target, $target_len, $target_start, $target_end, $strand, $scaffold_id, $start, $end, $score) = (split /\t/, $_)[0,1,2,3,4,5,7,8,10];
		my $coverage = sprintf "%.2f", 100*($target_end-$target_start+1)/$target_len;
		$target_end-$target_start+1 == $target_len or next;
		my @g = ();
		$g[0]{program} = 'BLASTN';
		$g[0]{type}    = 'gene';
		$g[0]{start}   = $start;
		$g[0]{end}     = $end;
		$g[0]{score}   = '.';
		$g[0]{strand}  = $strand;
		$g[0]{phase}   = '.';
		$g[0]{ann}     = "Target=$target $target_start $target_end;coverage=$coverage;";
		$g[1]{program} = 'BLASTN';
		$g[1]{type}    = 'CDS';
		$g[1]{start}   = $start;
		$g[1]{end}     = $end;
		$g[1]{score}   = $score;
		$g[1]{strand}  = $strand;
		$g[1]{phase}   = 0;
		$g[1]{ann}     = '';
		push @{$$gff{$scaffold_id}}, \@g;
	}
	close IPUT;
}


sub Read_Blast {
	my ($iput, $blast) = @_;
	open my $IPUT, "<$iput" or die;
	my @key;
	while (<$IPUT>) {
		/^# Fields: (.*)$/ or next;
		@key = split ', ', $1;
		last;
	}
	@key or die;
	my %title = ('query id'          => 'Query_ID',
	             'subject id'        => 'Subject_ID',
	             'query length'      => 'Query_len',
	             'subject length'    => 'Subject_len',
	             'evalue'            => 'Evalue',
	             '% identity'        => '%Ident',
	             '% positives'       => '%Pos',
	             'alignment length'  => 'Align_len',
	             'mismatches'        => 'Mismatch',
	             'gap opens'         => 'Gap_opening',
	             'q. start'          => 'Query_start',
	             'q. end'            => 'Query_end',
	             's. start'          => 'Subject_start',
	             's. end'            => 'Subject_end',
	             'bit score'         => 'Score',
	             'subject tax ids'   => 'Tax_ID',
	             'subject sci names' => 'Sci_Name',
	             'subject title'     => 'Subject_Title');
	for (my $i=0;$i<@key;$i++) {
		$key[$i] = $title{$key[$i]} || $key[$i];
	}
	while (<$IPUT>) {
		/^#/ and next;
		my $h = TOOLS::read_value(\@key, $_);
		push @$blast, $h;
	}
	close $IPUT;
}

############################################################

sub extract_cds {
	my $gff = shift;
	foreach my $g (@$gff) {
		my @g;
		$g[0] = $$g[0];
		foreach ( sort {$$a{start}<=>$$b{start} || $$a{end}<=>$$b{end}} @$g) {
			$$_{type} eq 'CDS' and push @g, $_;
		}
		@$g = @g;
	}
}


sub extract_fasta {
	my ($scaffold, $gff) = @_;
	foreach my $g (@$gff) {
		my %s;
		my ($len, $seq);
		if (@$g > 1) {
			for (my $i=1;$i<@$g;$i++) {
				$$g[$i]{type} eq 'CDS' or next;
				$len  = $$g[$i]{end} - $$g[$i]{start} + 1;
				$seq .= substr $$scaffold{seq}, $$g[$i]{start}-1, $len;
			}
		}
		elsif (@$g == 1) {
			$len = $$g[0]{end} - $$g[0]{start} + 1;
			$seq = substr $$scaffold{seq}, $$g[0]{start}-1, $len;
		}
		if ($$g[0]{strand} eq '-') {
			$seq = reverse $seq;
			$seq =~ tr/ACGTRYMK/TGCAYRKM/;
			$seq =~ tr/acgtrymk/tgcayrkm/;
		}
		$s{seq} = uc $seq;
		my $N = $seq =~ tr/Nn//;
		$s{gap} = $N;
		$s{ann} = "locus=$$scaffold{id}:$$g[0]{start}:$$g[0]{end}:$$g[0]{strand};";
		$s{ann}.= " $$g[0]{ann}" if ($$g[0]{ann});
		$$g[0]{fasta} = \%s;
	}
}


sub extract_fa {
	my ($gff, $scaffold) = @_;
	foreach my $g (@$gff) {
		my %s;
		my ($len, $seq);
		$len = $$g[0]{end} - $$g[0]{start} + 1;
		$seq = substr $$scaffold{seq}, $$g[0]{start}-1, $len;
		if ($$g[0]{strand} eq '-') {
			$seq = reverse $seq;
			$seq =~ tr/ACGTRYMK/TGCAYRKM/;
			$seq =~ tr/acgtrymk/tgcayrkm/;
		}
		$s{seq} = uc $seq;
		my $N = $seq =~ tr/Nn//;
		$s{gap} = $N;
		$s{ann} = "locus=$$scaffold{id}:$$g[0]{start}:$$g[0]{end}:$$g[0]{strand};";
		$s{ann}.= " $$g[0]{ann}" if ($$g[0]{ann});
		$$g[0]{fasta} = \%s;
	}
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


sub filter_gap {
	my ($gff, $N) = @_;
	$N ||= 0;
	foreach my $g (@$gff) {
		$$g[0]{filter} = 1 if ($$g[0]{fasta}{gap} > $N);
	}
}

__END__
