#adjust.pl:change the input species_abundance list and product the new-percentage.nctf file
use strict; use warnings;

my $old=$ARGV[0];
my $index=$ARGV[1];
my $line1;
my $line2;
my $flag;
my @array1;
my @name;
my @array2;
my $percentage;
my $total_percentage=0;
my $regular;
my $i;
my @array_record;
my $number=0;
my $usage="perl $0 input-composition-table index-file";
my $whole;

$i=@ARGV;
if($i!=2)
{
	print "$usage\n";
	exit;
}

open(OLD,$old) or die "Can't open the $old file!\nusage:$usage\n";
open(NEW,">new-percentage.txt");

while(<OLD>)
{
	chomp;
 	$line1=$_;
 	@array1=split("\t",$line1);
 	$percentage=$array1[1];

#new percentage=old percentage * length
 	open(INDEX,$index) or die "usage:$usage\n";
 	$flag=0;

	#change the pattern of the name
 	@name=split("",$array1[0]);
	$regular="";

 	for ($i=0;$i<@name;$i++)
 	{
  		if (($name[$i]!~/\w/)&&($name[$i]!~/\d/)&&($name[$i]!~/\s/))
  		{ 
   			$regular=$regular."\\".$name[$i];
  		}
  		else
  		{
   			$regular=$regular.$name[$i];
  		}
 	} 
	
	$whole=0;
 	while(<INDEX>)
 	{
  		chomp;
  		$line2=$_;
		@array2=split("\t",$line2);
  		if (($array2[3]=~/$regular/)&&($array2[2]==0))
  		{
   			$flag=1;
   			$whole+=$array2[1];
  		}
 	}
 	if ($flag==0)
 	{
  		print "In the index don't have the: $array1[0]!\n";
  		exit;
 	}
 	close INDEX;
 	$percentage*=$whole;
 	$total_percentage+=$percentage;
 	$array_record[$number][0]=$regular;
 	$array_record[$number][1]=$percentage;
	$array_record[$number][2]=$whole;
 	$number++;
}
close OLD;

#the final percentage after trimming
for ($i=0;$i<$number;$i++)
{
	open(INDEX,$index);
	while(<INDEX>)
	{
		chomp;
		$line2=$_;
		if($line2=~/$array_record[$i][0]/)
		{
			@array2=split("\t",$line2);
			if($array2[2]!=0)
			{
				next;
			}
			$percentage=$array_record[$i][1]/$total_percentage*$array2[1]/$array_record[$i][2];
			printf NEW ("%s\t%0.5f\n",$array2[3],$percentage);
		}
	}
	close INDEX;
}
close NEW;
