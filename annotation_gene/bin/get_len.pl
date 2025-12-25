open A,$ARGV[0] or die $!;
while(<A>)
{
	chomp;
	if(/>(\S+)/)
	{
		$id=$1;
	}
	else
	{
		$len{$id}+=length $_;
	}
}
close A;

foreach $id (sort {$len{$b}<=>$len{$a}} keys %len)
{
	print "$id\t$len{$id}\n";
}
