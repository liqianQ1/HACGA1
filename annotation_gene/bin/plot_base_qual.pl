#!/usr/bin/perl
use strict;
use warnings;
# use lib "/home/zhanghk/bins/qc";
use lib "/share/public/software/lib_ssinfo/perl/qc/";
use Term::ANSIColor qw(:constants);
   $Term::ANSIColor::AUTORESET=1;
use FindBin qw($Bin);
use Getopt::Long;
use PLOT;

my $usage=<<"USAGE";

	Script   : plot distribution of base and quality along reads from fq1.check and/or fq2.check file
	Usage    : perl $0 -1 <fastq1_file> -2 <fastq2_file> -o [output_file_prefix]
	Exmple   : perl $0 -1 example_fq1.check -2 example_fq2.check -o ./example

USAGE

print GREEN join ' ', ($0,@ARGV,"\n");

############################################################
my ($fq1, $fq2, $output_pre);
GetOptions(
	"1=s"=>\$fq1,
	"2=s"=>\$fq2,
	"o=s"=>\$output_pre,
);

die $usage unless ($fq1);
#my ($in, $output_pre) = @ARGV;

if (!$output_pre) {
	$output_pre = $fq1;
	$output_pre =~ s/check$//;
}
if ($output_pre !~ /[\/\.]$/) {
	$output_pre .= '.';
}
my $in = $fq1;
my $head = "type\tfq1";
my @qc;
if($fq2){
	$head.="\tfq2\tall";
	$in = $output_pre."temp";
	my $cycle = 0;
	open IN1,"<$fq1" or die "error $fq1\n";
	open OUT,">$in" or die "create error!\n";
	while(my $line = <IN1>){
		chomp $line;
		#print $line;
		if($line =~ /^\d/){
			print OUT "$line\n";
			$cycle++;
		}elsif($line !~ /^#/){
			push @qc,$line;
			print @qc."\n";
		}
	}
	close IN1;

	open IN2,"<$fq2" or die "error $fq2";
	my @arry = <IN2>;
	my @arry2 = reverse(@arry);
	my $index = 0;
	foreach my $line (@arry){
		chomp $line;
		if($line =~ /^\d+/){
			$cycle++;
			$line =~ s/^\d+/$cycle/;
			print OUT "$line\n";
		}elsif($line !~ /^#/){
			my ($t,$v) = split /\s+/,$line;
			$qc[$index].="\t$v";
			$index++;
		}
	}
	close IN2;
	close OUT;
}

print "output:\t${output_pre}base.png\n";
print "output:\t${output_pre}qual2.png\n";

my ($qual, $cycle, $qq) = read_fastqc ($in);
my $temp = $output_pre.".temp";
my $fqc =  $output_pre."cleandata.xls";
open QC,">$fqc";
my $oqc = @qc==0?$qq:\@qc;
print QC "$head\n";
foreach my $line (@{$oqc}){
	my @a = split /\t+/,$line;
	if(scalar(@a)>2){
		my $all = $a[1]+$a[2];
		if($a[0]=~/%/){
			$all=$all/2;
		}
		print QC "$a[0]\t$a[1]\t$a[2]\t$all\n";
	}else{
		print QC "$line\n";
	}
}
close QC;
my $out = "${output_pre}base";
my $plot;
$plot = "reset;\n".
        "set terminal postscript portrait color size 8, 5;\n".
        "set out      '$out.ps';\n".
        "set grid     front lc rgb \'gray\';\n".
        "set title    'Per base sequence content';\n".
        "set xlabel   'Position in read';\n".
        "set ylabel   'Percent';\n".
        "set xtics    10,20,$cycle;\n".
        "set xrange   [0:$cycle+1];\n".
        "set yrange   [-2:52];\n";
$plot.= "plot '$in' u 1:2 w l lw 3 t 'A'".
              ", '' u 1:3 w l lw 3 t 'C'".
              ", '' u 1:4 w l lw 3 t 'T'".
              ", '' u 1:5 w l lw 3 t 'G'".
              ", '' u 1:6 w l lw 3 t 'N'";
PLOT::plot ($plot, $out);

=cut
$out = "${output_pre}qual1";
$plot = "reset;\n".
        "set terminal postscript portrait color size 16, 5;\n".
	    "set out      '$out.ps';\n".
        "set grid     front lc rgb \'gray\';\n".
	    "set title    'Per base sequence quality';\n".
	    "set xlabel   'Position in read';\n".
		"set ylabel   'Quality';\n".
        "set xtics    10,10,$cycle;\n".
		"set xrange   [0:$cycle+1];\n".
		"set yrange   [-3:43];\n".
		"set ytics    0,10,40;\n".
		"set key      off;\n".
		"set boxwidth 0.7;\n";
$plot .= "plot '$in' u 1:10:12:9:13 w candlesticks whiskerbars lc rgb 'blue' lw 2".
               ", '' u 1:11:11:11:11 w candlesticks lc rgb 'red' lw 2".
               ", '' u 1:14 w l lc 'blue' lw 2;\n";
PLOT::plot ($plot, $out);


=cut


$out = "${output_pre}qual";
$plot = "reset;\n".
        "set terminal postscript portrait color size 8, 5;\n".
	    "set out      '$out.ps';\n".
        "set grid     front lc rgb \'gray\';\n".
		"set rmargin 10;\n".
	    "set title    'Per base sequence quality';\n".
	    "set xlabel   'Position in read';\n".
		"set ylabel   'Quality';\n".
        "set xtics    10,20,$cycle;\n".
		"set xrange   [0:$cycle+1];\n".
		"set yrange   [-3:43];\n".
		"set ytics    0,10,40;\n".
		"set key      off;\n".
		"set palette  defined (0 '\#ffffff', 10 '\#00ff00', 30 '\#ffff00', 50 '\#ff0000',  100 '\#800000');\n".
		"set font ',5';\n".
		"set style line 2604 linetype -1 linewidth .4;\n".
		"set colorbox noborder;\n".
		"set cbtics 20,20,100;\n".
		"set colorbox user origin 0.9,0.25 size .01,0.5;\n";
$plot .="plot '-' u 1:2:3 w image\n";
$plot .= (join "\t", @{$_})."\n" foreach (@{$qual});
$plot .= "e\n";
PLOT::plot ($plot, $out);
unlink($in);


sub read_fastqc {
	my $fastqc = shift;
	my @qual;
	my $cycle=0;
	my @qq;
	open FASTQC, "<$fastqc" or die "Error: cannot open $fastqc";
	while (my $line = <FASTQC>) {
		chomp $line;
		if ($line =~ /^\d+/) {
			$cycle++;
			my @col = split /\s+/, $line;
			for (my $i=6; $i<@col; $i++) {
				my $q = $i-6;
				my @tmp = ($cycle, $q, $col[$i]);
				push @qual, \@tmp;
			}
		}elsif($line !~ /^#/){
			push @qq, $line;
		}
	}
	close FASTQC;


	return \@qual, $cycle,\@qq;
}

__END__
