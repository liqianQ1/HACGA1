#!/usr/bin/perl
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
	Usage:	$0 <validate.log> <input> [<input.fa>] [-o <output>]
USAGE

print GREEN join ' ', ($0,@ARGV,"\n");

############################################################

my ($oput);
GetOptions ("o:s"=>\$oput);
die $usage unless (@ARGV>=2);
my ($log, $iput, $fa) = @ARGV;
$oput ||= 'out';


my $program;
open IPUT, "<$iput" or die;
while (<IPUT>) {
	/^LOCUS/ and $program='AUGUSTUS' and last;
	/^>/ and $program='SNAP' and last;
}
close IPUT;
$program or die;

my %badlist;
my %goodlist;
open LOG, "<$log" or die;
if ($program eq 'AUGUSTUS') {
	while(<LOG>) {
		/[Ss]equence (\S+):/ or next;
		$badlist{$1} = 1;
	}
}
elsif ($program eq 'SNAP') {
	while(<LOG>) {
		/(\S+)\s+OK/ or next;
		$goodlist{$1} = 1;
	}
}
close LOG;


my %scaflist;
open IPUT, "<$iput" or die;
open OPUT, ">$oput" or die;
if ($program eq 'AUGUSTUS') {
	$/="\n//\n";
	while(<IPUT>) {
		/^LOCUS +(\S+) .*/ or die;
		my $id = $1;
		print OPUT $_ if (!$badlist{$id});
	}
}
elsif ($program eq 'SNAP') {
	my $scafid;
	while(<IPUT>) {
		if (/^>(\S+)/) {
			$scafid = $1;
		}
		else {
			/(\S+)\n/ or die;
			if ($goodlist{$1}) {
				if ($scafid) {
					print OPUT ">$scafid\n";
					$scaflist{$scafid} = 1;
					$scafid = '';
				}
				print OPUT $_;
			}
		}

	}
}
close IPUT;
close OPUT;


if ($program eq 'SNAP') {
	$oput =~ s/\.\w+$//;
	$oput .= '.dna';
	open OPUT, ">$oput" or die;
	open IPUT, "<$fa";
	my $scafid;
	$/=">";
	<IPUT>;
	while (<IPUT>) {
#		$/="\n>";
		chomp;
		/^(\S+)/ and $scafid = $1;
		my @line = split/\n/,$_,2;
#		print "$1###$line[0]\n";
#		/^>(\S+)/ and $scafid = $1;
		if ($scaflist{$scafid}) {
			print OPUT ">$_";
		}
	}
	close IPUT;
	close OPUT;
	$/="\n";
}

__END__
