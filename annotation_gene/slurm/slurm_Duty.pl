#!/usr/bin/perl

=head1 Name

slurm-Duty.pl -- control processes running on linux slurm system

=head1 Description

This program throw the jobs and control them running on linux Slurm system. It reads jobs
from an input shell file. One line is the smallest unit of  a single job, however, you can also specify the
number of lines to form a single job. For sequential commands, you'd better put them
onto a single line, seperated by semicolon. In anywhere, "&" will be removed
automatically. The program will terminate when all its jobs are perfectly finished.

If you have so many jobs, the efficency depends on how many CPUs you can get,
or which queque you have chosen by --queue option. You can use the --maxjob option to
limit the number of throwing jobs, in order to leave some CPUs for other people.
When each job consumes long time, you can use the --interval option to increase interval
time for qstat checking , in order to reduce the burden of the head node.

As Slurm can only recognize absolute path, so you'd better use absolute path everywhere,
we have developed several ways to deal with path problems:
(1) We have added a function that converting local path to absolute
path automatically. If you like writting absolute path by yourself, then you'd better close this
function by setting "--convert no" option.
(2) Note that for local path, you'd better write
"./me.txt" instead of only "me.txt", because "/" is the  key mark to distinguish path with
other parameters.
(3) If an existed file "me.txt" is put in front of the redirect character ">",
or an un-created file "out.txt" after the redirect character ">",
the program will add a path "./" to the file automatically. This will avoid much
of the problems which caused by forgetting to write "./" before file name.
However, I still advise you to write "./me.txt" instead of just "me.txt", this is a good habit.
(4) Please also note that for the re-direct character ">" and "2>", there must be space characters
both at before and after, this is another good habit.

There are several mechanisms to make sure that all the jobs have been perfectly finished:
(1) We add an auto job completiton mark "This-Work-is-Completed!" to the end of the job, and check it after the job finished
(2) We check "GLIBCXX_3.4.9 not found" to make sure that the C/C++ libary on computing nodes are in good state
(3) We provide a "--secure" option to allow the users define their own job completition mark. You can print a mark
    (for example, "my job complete") to STDERR at the end of your program, and set --secure "my job complete" at
	this program. You'd better do this when you are not sure about wheter there is bug in your program.
(4) We provide a "--reslurm" option, to throw the unfinished jobs automatically, until all the jobs are
    really finished. By default, this option is closed, please set it forcely when needed. The maximum
	reslurm cycle number allowed is 1000.
(5) Add a function to detect the died computing nodes automatically.
(6) Add checking "iprscan: failed" for iprscan
(7) Add a function to detect queue status, only "r", "t", and "qw" is considered correct.
(8) Add check "failed receiving gdi request"

Normally, The result of this program contains 3 parts: (Note that the number 24137 is the process Id of this program)
(1) work.sh.24137.globle,     store the shell scripts which has been converted to global path
(2) work.sh.24137.slurm,       store the middle works, such as job script, job STOUT result, and job STDERR result
(3) work.sh.24137.log,      store the error job list, which has been throwed more than one times.

I advice you to always use the --reslurm option and check the .log file after this program is finished. If you find "All jobs finished!", then
then all the jobs have been completed. The other records are the job list failed in each throwing cycle, but
don't worry, they are also completed if you have used --reslurm option.

For the resource requirement, by default, the --resource option is set to vf=1.9G, which means the total
memory restriction of one job is 1.9G. By this way, you can throw 8 jobs in one computing node, because the
total memory restriction of one computing node is 15.5G. If your job exceeds the maximum memory allowed,
then it will be killed forcely. For large jobs, you must specify the --resource option manually, which
has the same format with "slurm -l" option. If you have many small jobs, and want them to run faster, you
also need to specify a smaller memory requirement, then more jobs will be run at the same time. The key
point is that, you should always consider the memory usage of your program, in order to improve the efficency
of the whole cluster.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Autor: Hu Yujie  huyj@genomics.org.cn
  Version: 8.2,  Date: 2009-11-11

=head1 Usage

  perl slurm-Duty.pl <jobs.txt>
  --global          only output the global shell, but do not excute
  --partition <str>     specify the partition to use, default all availabile queues
  --interval <num>  set interval time of checking by qstat, default 300 seconds
  --lines <num>     set number of lines to form a job, default 1
  --maxjob <num>    set the maximum number of jobs to throw out, default 100
  --convert <yes/no>   convert local path to absolute path, default yes
  --secure <mark>   set the user defined job completition mark, default no need
  --reslurm          reslurm the unfinished jobs untill they are finished, default no
  --mem <str>  set the required memory resource used in sbatch --mem option, default vf=1G
  --cpu <num>       set cpu num of per job,default 1
  --jobprefix <str> set the prefix tag for slurmed jobs, default work
  --verbose         output verbose information to screen
  --help            output help information to screen

  --getmem          output the usage (example: cpu=00:26:45, mem=111.63317 GBs, io=0.00000, vmem=259.148M, maxvmem=315.496M);

=head1 Exmple

  1.work with default options (the most simplest way)
  perl slurm-Duty.pl ./work.sh

  2.work with user specifed options: (to select queue, set checking interval time, set number of lines in each job, and set number of maxmimun running jobs)
  perl slurm-Duty.pl --partition node01 -interval 1 -lines 3 -maxjob 10  ./work.sh

  3.do not convert path because it is already absolute path (Note that errors may happen when convert local path to absolute path automatically)
  perl slurm-Duty.pl --convert no ./work.sh

  4.add user defined job completion mark (this can make sure that your program has executed to its last sentence)
  perl slurm-Duty.pl -inter 1  -secure "my job finish" ./work.sh

  5.reslurm the unfinished jobs until all jobs are really completed (the maximum allowed reslurm cycle is 10000)
  perl slurm-Duty.pl --reslurm ./work.sh

  6.work with user defined memory usage
  perl slurm-Duty.pl --mem 4G ./work.sh

  7.recommend combination of usages for common applications (I think this will suit for 99% of all your work)
  perl slurm-Duty.pl --partition all.q --mem 4G -maxjob 10 --reslurm ./work.sh

=cut

BEGIN{push(@INC,"/opt/conda/envs/ann/lib/perl5")};
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(rmtree);
use Cwd qw(abs_path);
use Data::Dumper;

##get options from command line into variables and set default values
my ($Global,$Queue, $Interval, $Lines, $Maxjob, $cpu, $Convert,$Secure,$Reslurm, $Reslurm2,$Resource,$Job_prefix,$Verbose, $Help, $getmem);
GetOptions(
	"global"=>\$Global,
	"lines:i"=>\$Lines,
	"maxjob:i"=>\$Maxjob,
    "cpu:i"=>\$cpu,
	"interval:i"=>\$Interval,
	"partition:s"=>\$Queue,
	"convert:s"=>\$Convert,
	"secure:s"=>\$Secure,
	"reslurm"=>\$Reslurm,
	"mem:s"=>\$Resource,
	"jobprefix:s"=>\$Job_prefix,
	"verbose"=>\$Verbose,
	"help"=>\$Help,
	"getmem"=>\$getmem,
);
$Queue ||= "q_cn";
$Interval = 120;
#$Interval = 300;
$Lines ||= 1;
$Maxjob ||= 100;
$cpu ||=1;
$Convert ||= 'yes';
$Resource ||= "1G";
$Job_prefix ||= "work";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $work_shell_file = shift;

##global variables
my $work_shell_file_globle = $work_shell_file.".$$.globle";
my $work_shell_file_error = $work_shell_file.".$$.log";
my $Work_dir = $work_shell_file.".$$.slurm";
my $current_dir = abs_path(".");

### get mem ;  add by nixiaoming nixiaoming@genomics.cn
my $work_shell_mem = $work_shell_file.".$$.mem.log";
my %meminfo=();
#open GETMEM,'>',$work_shell_mem or die "can't open the mem info $work_shell_mem ";

my $whoami=`whoami`;
chomp($whoami);
	#print GETMEM "User:\t\t$whoami\nShellPath:\t$current_dir/$Work_dir\n";
#close GETMEM;
###

if ($Convert =~ /y/i) {
	absolute_path($work_shell_file,$work_shell_file_globle);
}else{
	$work_shell_file_globle = $work_shell_file;
}

if (defined $Global) {
	exit();
}

my $time='`date +%F'."'  '".'%H:%M`'; ### add by nixiaoming

## read from input file, make the slurm shell files
my $line_mark = 0;
my $Job_mark="00001";
mkdir($Work_dir);
my @Shell;  ## store the file names of slurm sell
my $allline = `wc -l $work_shell_file_globle`;
chomp $allline;
open IN, $work_shell_file_globle || die "fail open $work_shell_file_globle";
while(<IN>){
	chomp;
	#s/&/;/g;
	next unless($_);
	if ($line_mark % $Lines == 0) {
		open OUT,">$Work_dir/$Job_prefix\_$Job_mark.sh" || die "failed creat $Job_prefix\_$Job_mark.sh";
        print OUT "#!/bin/bash\n";
		print OUT 'echo start at time '.$time."\n";
		push @Shell,"$Job_prefix\_$Job_mark.sh";
		$Job_mark++;
	}
	s/;\s*$//;  ##delete the last character ";", because two ";;" characters will cause error in slurm
	s/;\s*;/;/g;
    my $tmp=$_;
	if($_!~/&$/){
		print OUT $_.' &&  echo This-Work-is-Completed!'."\n";
	}else{
		chop $_;
		print OUT $_.' &&  echo '.'"'."$_".'"'.' This-Work-is-Completed! &'."\n";
	}
	#print OUT $_."; echo This-Work-is-Completed!\n";

	if (($line_mark % $Lines == $Lines - 1) || ($line_mark==$allline-1)) {
        if($tmp=~/&$/){
            print OUT "wait\n";
        }
		print OUT 'echo finish at time '.$time."\n";
		close OUT;
	}

	$line_mark++;
}
close IN;
close OUT;


print STDERR "make the slurm shell files done\n" if($Verbose);


## run jobs by slurm, until all the jobs are really finished
my $slurm_cycle = 1;
if(defined($getmem)){
	open GETMEM,'>',$work_shell_mem or die "can't open the mem info $work_shell_mem ";
	print GETMEM "User:\t\t$whoami\nShellPath:\t$current_dir/$Work_dir\n";
}
while (@Shell) {

	## throw jobs by slurm
	##we think the jobs on died nodes are unfinished jobs
	my %Alljob; ## store all the job IDs of this cycle
	my %Runjob; ## store the real running job IDs of this cycle
	my %Error;  ## store the unfinished jobs of this cycle
	chdir($Work_dir); ##enter into the slurm working directoy
	my $job_cmd = "sbatch --mem $Resource -c $cpu ";  ## -l h_vmem=16G,s_core=8
	$job_cmd .= "-p $Queue "  if(defined $Queue); ##set queue
	die "no valid resource\n" unless($Resource);
	my $part_resource = $Resource;
	#$job_cmd .= "-l $Resource " if(defined $Resource); ##set resource
	##warn $job_cmd;
	my $job_id;
	for (my $i=0; $i<@Shell; $i++) {
		while (1) {
			my $run_num = run_count(\%Alljob,\%Runjob,\%meminfo);
			if ( $run_num != -1 && ($i < $Maxjob || $run_num < $Maxjob) ) {#此处修复了一处问题：如果run_num==-1且i < Maxjob，集群有问题情况下还是会投递任务，那么后面将出现未知情况
				my $jod_return = `$job_cmd -o $Shell[$i].o%j -e $Shell[$i].e%j $Shell[$i]`;
				sleep(2);#20160219 投递会漏任务，可能是间隔太短，因此休眠2秒再投
				if($jod_return =~ /Submitted batch job (\d+)/){
					$job_id = $1;
				}else{
					warn "can't get job id from: $jod_return\n";
				}
				$Alljob{$job_id} = $Shell[$i];  ## job id => shell file name
				print STDERR "throw job $job_id in the $slurm_cycle cycle\n" if($Verbose);	
				last;
			}else{
				print STDERR "wait for throwing next job in the $slurm_cycle cycle\n" if($Verbose);
				sleep $Interval;
			}
		}
	}
	chdir($current_dir); ##return into original directory


	###waiting for all jobs fininshed
	while (1) {
		my $run_num = run_count(\%Alljob,\%Runjob,\%meminfo);
		last if($run_num == 0);
		print STDERR "There left $run_num jobs runing in the $slurm_cycle cycle\n" if(defined $Verbose);
		
		if($run_num != -1 and defined $getmem){ ### get mem ;  add by nixiaoming nixiaoming@genomics.cn
			foreach my $shname (sort keys %meminfo){
				my $jobinfo=$meminfo{$shname};
				chomp $jobinfo;
				$jobinfo =~ s/usage\s*\w*:\s*//g;
				$time=`date`;
				print GETMEM "$time\t$whoami\t$shname\t$jobinfo\n";
			}
		}
		
		sleep $Interval;
	}

	print STDERR "All jobs finished, in the firt cycle in the $slurm_cycle cycle\n" if($Verbose);


	##run the secure mechanism to make sure all the jobs are really completed
	open OUT, ">>$work_shell_file_error" || die "fail create $$work_shell_file_error";
	chdir($Work_dir); ##enter into the slurm working directoy
	foreach my $job_id (sort keys %Alljob) {
		my $shell_file = $Alljob{$job_id};

		##read the .o file
		my $content;
		if (-f "$shell_file.o$job_id") {
			#open IN,"$shell_file.o$job_id" || warn "fail $shell_file.o$job_id";
			#$content = join("",<IN>);
			#close IN;
			$content = `tail -n 1000 $shell_file.o$job_id`;
		}
		##check whether the job has been killed during running time
		if ($content !~ /This-Work-is-Completed!/) {
			$Error{$job_id} = $shell_file;
			print OUT "In slurm cycle $slurm_cycle, In $shell_file.o$job_id,  \"This-Work-is-Completed!\" is not found, so this work may be unfinished\n";
		}


		##read the .e file
		my $content;
		if (-f "$shell_file.e$job_id") {
			#open IN,"$shell_file.e$job_id" || warn "fail $shell_file.e$job_id";
			#$content = join("",<IN>);
			#close IN;
			$content = `tail  -n 1000 $shell_file.e$job_id`;
		}
		##check whether the C/C++ libary is in good state
        if ($content =~ /oom-kill/){
            $Error{$job_id} = $shell_file;
        }
		if ($content =~ /GLIBCXX_3.4.9/ && $content =~ /not found/) {
			$Error{$job_id} = $shell_file;
			print OUT "In slurm cycle $slurm_cycle, In $shell_file.e$job_id,  GLIBCXX_3.4.9 not found, so this work may be unfinished\n";
		}

		##check whether iprscan is in good state
		if ($content =~ /iprscan: failed/) {
			$Error{$job_id} = $shell_file;
			print OUT "In slurm cycle $slurm_cycle, In $shell_file.e$job_id, iprscan: failed , so this work may be unfinished\n";
		}

		##check the user defined job completion mark
		if (defined $Secure && $content !~ /$Secure/) {
			$Error{$job_id} = $shell_file;
			print OUT "In slurm cycle $slurm_cycle, In $shell_file.o$job_id,  \"$Secure\" is not found, so this work may be unfinished\n";
		}
	}

	##make @shell for next cycle, which contains unfinished tasks
	@Shell = ();
	foreach my $job_id (sort keys %Error) {
		my $shell_file = $Error{$job_id};
		push @Shell,$shell_file;
	}

	$slurm_cycle++;
	if($slurm_cycle > 10){
		print OUT "\n\nProgram stopped because the reslurm cycle number has reached 10, the following jobs unfinished:\n";
		foreach my $job_id (sort keys %Error) {
			my $shell_file = $Error{$job_id};
			print OUT $shell_file."\n";
		}
		print OUT "Please check carefully for what errors happen, and redo the work, good luck!";
		die "\nProgram stopped because the reslurm cycle number has reached 10000\n";
	}

	print OUT "All jobs finished!\n" unless(@Shell);

	chdir($current_dir); ##return into original directory
	close OUT;
	print STDERR "The secure mechanism is performed in the $slurm_cycle cycle\n" if($Verbose);

	last unless(defined $Reslurm);
}

if(defined $getmem){ ### get mem ;  add by nixiaoming nixiaoming@genomics.cn
	print GETMEM "User:\t\t$whoami\nShellPath:\t$current_dir/$Work_dir\n";
	foreach my $shname (sort keys %meminfo){
		my $jobinfo=$meminfo{$shname};
		chomp $jobinfo;
		$jobinfo =~ s/usage\s*\w*:\s*//g;
		print GETMEM "$whoami\t$shname\t$jobinfo\n";
	}
}

close GETMEM;
print STDERR "\nslurm-Duty.pl finished\n" if($Verbose);


####################################################
################### Sub Routines ###################
####################################################

sub absolute_path{
	my($in_file,$out_file)=@_;
	my($current_path,$shell_absolute_path);

	#get the current path ;
	$current_path=abs_path(".");

	#get the absolute path of the input shell file;
	if ($in_file=~/([^\/]+)$/) {
		my $shell_local_path=$`;
		if ($in_file=~/^\//) {
			$shell_absolute_path = $shell_local_path;
		}
		else{$shell_absolute_path="$current_path"."/"."$shell_local_path";}
	}

	#change all the local path of programs in the input shell file;
	open (IN,"$in_file");
	open (OUT,">$out_file");
	while (<IN>) {
		chomp;
		##s/>/> /; ##convert ">out.txt" to "> out.txt"
		##s/2>/2> /; ##convert "2>out.txt" to "2> out.txt"
		my @words=split /\s+/, $_;

		##improve the command, add "./" automatically
		for (my $i=1; $i<@words; $i++) {
			if ($words[$i] !~ /\//) {
				if (-f $words[$i]) {
					$words[$i] = "./$words[$i]";
				}elsif($words[$i-1] eq ">" || $words[$i-1] eq "2>"){
					$words[$i] = "./$words[$i]";
				}
			}

		}
		for (my $i=0;$i<@words ;$i++) {
			if (($words[$i]!~/^\//) && ($words[$i]=~/\//)) {
				$words[$i]= "$shell_absolute_path"."$words[$i]";
			}
		}
		print OUT join("  ", @words), "\n";
	}
	close IN;
	close OUT;
}


##get the IDs and count the number of running jobs
##the All job list and user id are used to make sure that the job id belongs to this program
##add a function to detect jobs on the died computing nodes.
sub run_count {
	my $all_p = shift;
	my $run_p = shift;
	my $memlist = shift;
	my $run_num = 0;

	%$run_p = ();
	#my $user = $ENV{"USER"} || $ENV{"USERNAME"};
	#my $user = $whoami;
	my $qstat_result = `squeue 2>\&1`;
	#$user = substr($user, 0, 12);
	if ($qstat_result =~ /failed receiving gdi request/ || $qstat_result =~ /error/) {
        sleep 120;
		$run_num = -1;
		return $run_num; ##系统无反应
	}
	my @jobs = split /\n/,$qstat_result;
	my %died;
	died_nodes(\%died) if (@jobs > 0); ##the compute node is down, 有的时候节点已死，但仍然是正常状态
	foreach my $job_line (@jobs) {
		$job_line =~s/^\s+//;
		my @job_field = split /\s+/,$job_line;
#		next if($job_field[3] ne $user);
		if (exists $all_p->{$job_field[0]}){
			my $node_name = $1 if($job_field[7] =~ /(comput\d+)/);
            #print($node_name);
            #print($job_field[4]);
            #if ( !exists $died{$node_name} && ($job_field[4] eq "PD" || $job_field[4] eq "R" || $job_field[4] eq "t" || $job_field[4] eq "s" || $job_field[4] eq "S" ) ) {
            if ( $job_field[4] eq "PD" || $job_field[4] eq "R" || $job_field[4] eq "t" || $job_field[4] eq "s" || $job_field[4] eq "S") {
				$run_p->{$job_field[0]} = $job_field[2]; ##job id => shell file name
				$run_num++;
				if ((defined $getmem) && ($job_field[4] eq "R")){### get mem ;  add by nixiaoming nixiaoming@genomics.cn
					my $jobinfo=`sstat $job_field[0] 2>&1 |sed -n '3p' `;
					$$memlist{$all_p->{$job_field[0]}}=$jobinfo;
				}
			}else{
				print STDERR "scancel $job_field[0] for job run in died node $node_name\n";#添加了打印语句 20160830
                #`scancel $job_field[0]`;
			}
		}
	}

	return $run_num; ##qstat结果中的处于正常运行状态的任务，不包含那些在已死掉节点上的僵尸任务
}


##HOSTNAME                ARCH         NCPU  LOAD  MEMTOT  MEMUSE  SWAPTO  SWAPUS
##compute-0-24 lx26-amd64 8 - 15.6G - 996.2M -
sub died_nodes{
	my $died_p = shift;

	my @lines = split /\n/,`sinfo`;
	shift @lines; ##remove the first three title lines

	foreach  (@lines) {
		my @t = split /\s+/;
		my $node_name = $t[0];
		my $memory_use = $t[5];
		$died_p->{$node_name} = 1 if($t[3]=~/-/);
	}

}
