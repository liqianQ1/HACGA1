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
	Usage:	$0 <scaffold.masked.fa> <glimmerhmm.gff> [-o <output.gff>]
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
	ANNOTATION::Read_GFF ($_, \%gff);
}


my $count;
open my $OPUT, ">$oput_gff" or die;
open my $IPUT, "<$iput" or die;
while (!eof $IPUT) {
	my $scaffold = FASTA::read_fasta_one_seq ($IPUT);
	my $gff = $gff{$$scaffold{id}};
	$gff or next;
	motify ($_) foreach (@$gff);
	ANNOTATION::sort_gff ($gff);
	ANNOTATION::extract_fa ($gff, $scaffold);
#	ANNOTATION::filter_gap ($gff, 1000);
	ANNOTATION::unique_gene_id ($gff, \$count, $$scaffold{id}, 'Genewise');
	ANNOTATION::Print_GFF ($OPUT, $gff, $$scaffold{id});
}
close $OPUT;
close $IPUT;


sub motify {
    my $g = shift;
    my @g;
    
    # 跳过 gene 类型的行
    return if $$g[0]{type} eq 'gene';
    
    # 如果是 mRNA 类型，则处理为 gene 类型
    if ($$g[0]{type} eq 'mRNA') {
        $g[0]{program} = $$g[0]{program};
        $g[0]{type}    = 'gene';
        $g[0]{start}   = $$g[0]{start};
        $g[0]{end}     = $$g[0]{end};
        $g[0]{score}   = '.';
        $g[0]{strand}  = $$g[0]{strand};
        $g[0]{phase}   = '.';
    }
    
    push @g, $$g[0];
    
    foreach (my $i = 1; $i < @$g; $i++) {
        # 跳过 gene 类型的记录
        next if $$g[$i]{type} eq 'gene';
        
        # 如果是 CDS 类型且前一个不是 exon 类型，则新增 exon 记录
        if ($$g[$i]{type} eq 'CDS' && $$g[$i-1]{type} ne 'exon') {
            my $j = @g;
            $g[$j]{program} = $$g[$i]{program};
            $g[$j]{type}    = 'exon';
            $g[$j]{start}   = $$g[$i]{start};
            $g[$j]{end}     = $$g[$i]{end};
            $g[$j]{score}   = '.';
            $g[$j]{strand}  = $$g[$i]{strand};
            $g[$j]{phase}   = '.';
        }
        push @g, $$g[$i];
    }

    # 清除注释字段
    $$_{ann} = '' foreach (@g);
    
    # 更新原始 GFF 数组
    @$g = @g;
}

__END__
