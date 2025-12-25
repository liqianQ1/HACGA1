#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage=<<USAGE;
	$0 -g 1.gff,2.gff, ... -o stat.out
USAGE

print GREEN join ' ', ($0,@ARGV,"\n");
###################################################################
my ($gff_list,$out);
GetOptions (
        "g:s"=>\$gff_list,
);
my @gff = split /\,/,$gff_list;

my (%stat);
foreach (@gff){
	my ($sam) = `basename $_`=~/([^\.]*)\./;
	stat_gff($_,$sam,\%stat);
}
print "\tNumber of genes\tMean CDS length(bp)\tExons per transcript\tMean exon length(bp)\tMean intron length(bp)\n";
foreach my $id(keys %stat){


	
	print "$id\t";
	printf "%d\t",$stat{$id}{gene}{num};
	printf "%d\t",$stat{$id}{CDS}{len}/$stat{$id}{mRNA}{num};
	printf "%.1f\t",$stat{$id}{exon}{num}/$stat{$id}{mRNA}{num};
	printf "%d\t",$stat{$id}{exon}{len}/$stat{$id}{exon}{num};
	printf "%d\n",$stat{$id}{intron}{len}/$stat{$id}{intron}{num};
}


##################################################################
sub stat_gff {
	my ($gff,$id,$hash)=@_;
	open GFF,"<$gff" or die $!;
	while (<GFF>){
		chomp ;
		next if ($_=~/^\s*&/);
 		my @line = split /\t/,$_;
		my $type = $line[2];
		my $len = $line[4]-$line[3]+1;
		${$hash}{$id}{$type}{num}++;
		${$hash}{$id}{$type}{len}+=$len;
	}
	close GFF;
	$$hash{$id}{intron}{num} = $stat{$id}{exon}{num}-$stat{$id}{mRNA}{num};
	$$hash{$id}{intron}{len} = $stat{$id}{mRNA}{len}-$stat{$id}{exon}{len};
} 


