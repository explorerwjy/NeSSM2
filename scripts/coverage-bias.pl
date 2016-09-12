####usage:perl $0 index composition-table bwa-file
use strict; use warnings;

my $line;
my @array;
my $gi;
my $my_gi;
my @result;
my $i;
my $pos;
my @a;
my $need;
my $line1;  ###composition table
my @array1;
my $line2;    ###index 
my @array2;
my $max;
my $total=0;
my $whole=0;
my $percentage;
my $flag;
my $usage="perl $0 index-file  composition-table  bwa-result-file\n";

#####check input
if(@ARGV!=3)
{
	print "ERROR!\n$usage";
	exit;
}

open(OUT,">coverage.txt") or die "Can't create the coverage.txt file\n";
open(IN,"<$ARGV[1]") or die "ERROR: Can't open the $ARGV[1] file!\n$usage";
while(<IN>)
{
	$line1=$_;
	@array1=split("\t",$line1);
	
	open(INDEX,"<$ARGV[0]") or die "ERROR: Can't open the $ARGV[0] file!\n$usage";
	$flag=0;
	while(<INDEX>)
	{
		$line2=$_;
		@array2=split("\t",$line2);
		if(($array2[3] eq $array1[0])&&($array2[2]==0))
		{
			$max=($array2[1]-($array2[1]%100))/100;
			$my_gi=$array2[0];
			$flag=1;
			last;
		}
	}
	close INDEX;
	if($flag==0)
	{
		print "Don't have the length information about $array1[0]!\n";
		exit;
	}

	for($i=0;$i<=$max;$i++)
	{
		$result[$i]=1;
	}

	open(BWA,"<$ARGV[2]") or die "ERROR: Can't open the $ARGV[2] file!\n$usage";
	while(<BWA>)
	{
		$line=$_;
		@array=split("\t",$line);
	
		if(@array<10)
		{
			next;
		}

		if($array[2] eq "*")
		{
			next;
		}

		($gi)=$array[2]=~/gi\|(\d+)\|/;
		if($gi==$my_gi)
		{
			$pos=($array[3]-($array[3]%100))/100;
			$result[$pos]++;
		}

		if($line!~/XA\:Z\:/)
		{
			next;
		}

		($need)=$line=~/XA\:Z\:(.+)\;$/;
		@a=split(";",$need);
		for($i=0;$i<@a;$i++)
		{
			($gi)=$a[$i]=~/gi\|(\d+)\|/;
			if($gi!=$my_gi)
			{
				next;
			}

			($pos)=$a[$i]=~/\|\,.(\d+)\,/;
			$pos=($pos-($pos%100))/100;
			$result[$pos]++;
		}
	}
	close BWA;

#####handle data
	$total=0;
	for($i=0;$i<=$max;$i++)
        {
                $total+=$result[$i];
        }
	print OUT "$my_gi=";
	$whole=$result[0];
	$percentage=$whole/$total;
	printf OUT ("%0.8f",$percentage);

	for($i=1;$i<=$max;$i++)
	{
		$whole+=$result[$i];
		$percentage=$whole/$total;
		printf OUT (":%0.8f",$percentage);
	}
	print OUT "\n";	
}
close IN;
close OUT;
