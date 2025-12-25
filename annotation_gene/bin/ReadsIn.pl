#!/usr/bin/perl -w
# Program Date:   2016.10.25

use strict;
use File::Basename qw(basename dirname);
use Getopt::Long;

my ($indir,$outdir);
GetOptions(
          "help|?" =>\&USAGE,
          "i:s"=>\$indir,
          "o:s"=>\$outdir,
);

if (!$indir ){
	print "perl $0 -i <inputdir>\n";
	exit;
}
$outdir||=$indir;
my @fq1=glob "$indir/*cleandata.xls";
open OUT,">$outdir/ReadsQuality.xls";

print OUT "Sample Name\tClean Reads\tClean Bases\tCleanN(%)\tCleanGC(%)\tCleanQ20(%)\tCleanQ30(%)\n";


foreach my $file(@fq1){
	my $inputdir= dirname $file;
	my $filename=basename $file;
	$filename=~s/\.cleandata\.xls//;
	print OUT "$filename";

	my ($readsnum1,$readsnum2,$basenum1,$basenum2,$Q201,$Q202,$Q301,$Q302,$N1,$N2,$GC1,$GC2)=(0,0,0,0,0,0,0,0,0,0,0,0);
	open QUAL,"$file";
	<QUAL>;
	while (<QUAL>){
		chomp;
		my @tem=split /\t/;
		if ($tem[0] ne "bases_number" && $tem[0] ne "nbases_number" && $tem[0] !~ /^AT/ &&  $tem[0] ne "reads_number"){
			if (@tem>3){
				print OUT "\t$tem[1];$tem[2]";
			}else{
				print OUT "\t$tem[1]";
			}
		}elsif($tem[0] eq "reads_number"){
			print OUT "\t$tem[1]";
		}elsif($tem[0] eq "bases_number"){
			$basenum1=$tem[1];
			if (@tem>3){
				my $base=$tem[1]*2;
				print OUT "\t$base";
				$basenum2=$tem[2];
			}else{
				print OUT "\t$tem[1]";
			}
		}elsif($tem[0] eq "nbases_number"){
			if (@tem>3){
				printf OUT "\t%4.2f;%4.2f",$tem[1]/$basenum1*100,$tem[2]/$basenum2*100;
			}else{
				printf OUT "\t%4.2f",$tem[2]/$basenum2*100;
			}
		}
	}
	print OUT "\n";
	close QUAL;
}

