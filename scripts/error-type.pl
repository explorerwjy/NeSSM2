####################################################
##  This program can get errors type and their 
##   probability from the .sam file. The file can be
##   created by aligning the fastq file to the database
##   through the BWA soft.
##
##   usage:perl $0 sam_file 454-or-illumina
####################################################
use strict; use warnings;

my $line;
my @array;
my @total;   #$total[0]:insertions,$total[1]:deletions,$total[2]:substitutions
my @temp;    
my @sequence;   #the read
my @result;   #reserve the substitution informations
my $i;
my $j;
my $insertion;   #how many insertions in one read
my $deletion;
my $substitution;
my $error_number;  #how many errors in one read
my $number;
my $get;  #reserve some results that come from the regular expression
my $position;
my $regular;
my $regular_get;  #get some value from regular expression
my $md;   #the MD tag
my $whole;
my $percentage;
my $usage="perl $0 sam_file 454-or-illumina";
my $need;   ##whether "N" in MD
my $a;
my $b;

$j=@ARGV;
if($j!=2)
{
	print "error!\n$usage\n";
	exit;
}

for($i=0;$i<3;$i++)    #initialize
{
	$total[$i]=0;
}

for($i=0;$i<4;$i++)
{
	for($j=0;$j<4;$j++)
	{
		$result[$i][$j]=0;
	}
}

if(($ARGV[1] ne "454")&&($ARGV[1] ne "illumina"))
{
	print "Please input the platform type: 454 or illumina\n";
	print "error!\n$usage\n";
	exit;
}

open(IN,$ARGV[0]) or die "can't open the $ARGV[0] file\n$usage\n";
while(<IN>)
{
	chomp;
	$line=$_;
	if($line=~/different NM/)   ##cause by samtools
        {
                next;
        }
	if($line=~/XA\:Z\:/)
	{
		($line)=$line=~/^(.+)\tXA\:Z\:/;
	}   ###use the first MD
	@array=split("\t",$line);
	if(@array>10)   #delete header
	{
		if($array[2] ne "*")       #delete the unmatched reads    #if-match read
		{
			($error_number)=$line=~/NM\:i\:(\d+)/;
			if($error_number==0)
			{
				next;
			}
			###if have "N" or other not "ATCG" in MD, exclude this line
			($need)=$line=~/MD\:Z\:(.+)$/;
			$a=$need=~tr/[A-Z]/[A-Z]/;
			$b=$need=~tr/ATCG/ATCG/;
			if($a!=$b)
			{
				next;
			} 

			$number=$array[5]=~tr/I/I/;    #how many times the "I" take place
			$insertion=0;
			
			if($number!=0)   #in this reads,there are insertions
			{
				@temp=split("I",$array[5]);
				for($i=0;$i<$number;$i++)
				{
					if($temp[$i]=~/[SMD]/)
					{
						($get)=$temp[$i]=~/[SMD](\d+)$/;
						$insertion+=$get;
					}
					else
					{
						$insertion+=$temp[$i];
					}
				}
			}
			$total[0]+=$insertion;

			$number=$array[5]=~tr/D/D/;     #in this read how many deletions
		 	$deletion=0;
			
			if($number!=0)
			{
                                @temp=split("D",$array[5]);
                                for($i=0;$i<$number;$i++)
                                {
                                        if($temp[$i]=~/[SMI]/)
                                        {
                                                ($get)=$temp[$i]=~/[SMI](\d+)$/;
                                                $deletion+=$get;
                                        }
                                        else
                                        {
                                                $deletion+=$temp[$i];
                                        }
                                }
                        }
                        $total[1]+=$deletion;		

			$substitution=$error_number-$insertion-$deletion;    #in this read, how many substitutions
			$total[2]+=$substitution;
						
			###########
			if($substitution>0)    #get the substitution informations     #if-1
			{
				@sequence=split("",$array[9]);
				if($array[5]=~/^\d+S/)     #when the soft clipping occur in the header of the read
				{
					($get)=$array[5]=~/^(\d+)S/;
					splice(@sequence,0,$get);    #delete the soft clipping
					($array[5])=$array[5]=~/^\d+S(.+)$/;
				}

				if($insertion>0)    #delete the insertions     #if-2
				{
					$position=0;
					$number=$array[5]=~tr/I/I/;
					$regular="";
					$i=0;
					$position=0;

					while($i!=$number)
					{
						if($array[5]=~/^$regular\d+M/)
						{
							$regular_get=$regular."(\\d+)M";
							($get)=$array[5]=~/^$regular_get/;
							$regular=$regular."\\d+M";
							$position+=$get;
						}
						elsif($array[5]=~/^$regular\d+D/)
						{
							$regular=$regular."\\d+D";
						}
						else
						{
							$regular_get=$regular."(\\d+)I";
							($get)=$array[5]=~/^$regular_get/;
							$i++;
							$regular=$regular."\\d+I";
							splice(@sequence,$position,$get);
						}
					}
				}           #end if-2
				
				($md)=$line=~/MD\:Z\:([\[0-9\]ATCG\^]+)/;
				if($deletion>0)    #delete deletion information in MD tag   #if-3
				{
					$md=~tr/^/Q/;
					@temp=split("Q",$md);
					$md=$temp[0];
					for($i=1;$i<@temp;$i++)
					{
						if($md=~/[ATCG]/)        #get the length before the deletion
						{
							($position)=$md=~/[ATCG](\d+)$/;
							($md)=$md=~/^(.*[ATCG])\d+$/;
						}
						else
						{
							$position=$md;
							$md="";
						}
					
						if($temp[$i]=~/\d+[ATCG]/)    #get the length after the deletion
						{
							($get)=$temp[$i]=~/^[ATCG]+(\d+)[ATCG]/;
							$position+=$get;
							($get)=$temp[$i]=~/^[ATCG]+\d+(.*)$/;
							$md=$md.$position.$get;
						}
						else
						{
							($get)=$temp[$i]=~/^[ATCG]+(\d+)$/;
							$position+=$get;
							$md=$md.$position;
						}
					}
				}   #end if-3

				##which base is changed by what
				$position=0;
				$number=$md=~tr/ATCG/ATCG/;
				if($number!=$substitution)    #test
				{
				##	print "This is an error alignment:$line\n";
					next;
				}
				
				$number=0;
				$regular="";
				while($number!=$substitution)
				{
					$regular_get=$regular."(\\d+)[ATCG]";
					($get)=$md=~/^$regular_get/;
					$position+=$get;
					$regular_get=$regular."\\d+(\\w)";
					($get)=$md=~/^$regular_get/;   #the orignal base
					$regular=$regular."\\d+[ATCG]";
					$number++;
					
					if($sequence[$position] eq "N")   #find the $i and $j in the @result
					{
						$total[2]--;
						$j=4;
					}
					elsif($sequence[$position] eq "A")
					{
						$j=0;
					}
                                        elsif($sequence[$position] eq "T")
                                        {
                                                $j=1;
                                        }
                                        elsif($sequence[$position] eq "C")
                                        {
                                                $j=2;
                                        }
                                        elsif($sequence[$position] eq "G")
                                        {
                                                $j=3;
                                        }
					$position++;

					if($get eq "A")
					{
						$i=0;
					}
					elsif($get eq "T")
					{
						$i=1;
					}
                                        elsif($get eq "C")
                                        {
                                                $i=2;
                                        }
                                        elsif($get eq "G")
                                        {
                                                $i=3;
                                        }
					if($j!=4)
					{
						$result[$i][$j]++;
					}
				}
			}    #end if-1
		}   #end if-match read		
	}
}
close IN;

open(OUT,">self-simulation.config") or die "can't create the self-simulation.config file\n";
print OUT "#substitution:(substitution+insertion):(substitution+insertion+deletion)\n";
print OUT $ARGV[1],"_type_probability_self=";
$whole=$total[0]+$total[1]+$total[2];
$percentage=$total[2]/$whole;
printf OUT ("%0.5f:",$percentage);
$percentage=($total[2]+$total[0])/$whole;
printf OUT ("%0.5f:1\n",$percentage);

$whole=0;
for($i=0;$i<4;$i++)
{
	for($j=0;$j<4;$j++)
	{
		$whole+=$result[$i][$j];
	}
}

print OUT "#AC:AG:AT:CG:CT:CA:GT:GA:GC:TA:TC:TG\n";
print OUT $ARGV[1],"_sub_probability_self=";
$percentage=$result[0][2]/$whole;
printf OUT ("%0.5f:",$percentage);

$percentage=$result[0][3]/$whole;
printf OUT ("%0.5f:",$percentage);

$percentage=$result[0][1]/$whole;
printf OUT ("%0.5f:",$percentage);

$percentage=$result[2][3]/$whole;
printf OUT ("%0.5f:",$percentage);

$percentage=$result[2][1]/$whole;
printf OUT ("%0.5f:",$percentage);

$percentage=$result[2][0]/$whole;
printf OUT ("%0.5f:",$percentage);

$percentage=$result[3][1]/$whole;
printf OUT ("%0.5f:",$percentage);

$percentage=$result[3][0]/$whole;
printf OUT ("%0.5f:",$percentage);

$percentage=$result[3][2]/$whole;
printf OUT ("%0.5f:",$percentage);

$percentage=$result[1][0]/$whole;
printf OUT ("%0.5f:",$percentage);

$percentage=$result[1][2]/$whole;
printf OUT ("%0.5f:",$percentage);

$percentage=$result[1][3]/$whole;
printf OUT ("%0.5f\n",$percentage);

#add informations of quality and read-length 
open(IN,"<quality.txt") or die "can't open the quality.txt\n";
while(<IN>)
{
	$line=$_;
	print OUT $line;
}
close IN;
close OUT;

