#usage:perl $0 fastq-file 454-or-illumina-or-pacbio
use strict; use warnings;

my $step;  #the quality range
my %quality;
my %ascii;
my $i;
my $j;
my $line;
my @array;
my $max=0;  #the max length of the reads
my $line_number=0;
my $length;
my $flag;
my $whole=0;
my $average;  #length average
my $sd;
my $number=0;  #how many reads accord with the length <=$max
my @part2;    #the second part's result
my $array_length;  #the part2's length
my $begin;
my %rand_length;
my $method;
my $usage="perl $0 fastq-file 454-or-illumina-or-pacbio";

$step=@ARGV;
if($step!=2)
{
	print "error!\n$usage\n";
	exit;
}
$step=72;
if($ARGV[1] eq "illumina")
{
	$method=1;
}
elsif($ARGV[1] eq "454")
{
	$method=2;
}
elsif($ARGV[1] eq "pacbio")
{
	$method=3;
}
else
{
	print "error!\n$usage\n";
	exit;
}

###############read the ascii.txt file
open(ASCII,"<ascii.txt") or die "Can't open the ascii.txt file!\n";
while(<ASCII>)
{
	chomp;
	$line=$_;
	@array=split("\t",$line);
	$ascii{$array[0]}=$array[1];
}
close ASCII;

open(IN,"<$ARGV[0]") or die "Can't open the $ARGV[0] file!\n$usage\n";
while(<IN>)
{
	chomp;
	$line_number++;
	$line=$_;

	if(($line_number%4)!=0)
	{
		next;
	}

	@array=split("",$line);
	$length=@array;
	if($length>$max)
	{
		$max=$length;
	}

	for($i=0;$i<$length;$i++)
	{
		statistics($array[$i],$i);
	}
}
close IN;

for($i=0;$i<$step*$max;$i++)  #absent add 0
{
	if(exists $quality{$i})
	{
		next;;
	}
	else
	{
		$quality{$i}=0;
	}
}

for($i=0;$i<$max;$i++)  #find the max length accord with the constrains and adjust quality
{
	$flag=adjust($i);
	if($flag==0)
	{
		$max=$i;   ##the max length
		last;
	}
}

$line_number=0;
open(IN,"<$ARGV[0]");
while(<IN>)   ###the average length
{
	chomp;
	$line=$_;
	if(($line_number%4)==3)
	{
		$length=length($line);
		if($length<=$max)
		{
			$whole+=$length;
			$number++;
		}
	}
	$line_number++;
}
close IN;
$average=$whole/$number;

$whole=0;
$line_number=0;
open(IN,"<$ARGV[0]");
while(<IN>)   ###the length's sd
{
        chomp;
        $line=$_;
        if(($line_number%4)==3)
        {
                $length=length($line);
                if($length<=$max)
                {
                        $whole+=($length-$average)**2;
			if(exists $rand_length{$length})
			{
				$rand_length{$length}++;
			}
			else
			{
				$rand_length{$length}=1;
			}
                }
        }
        $line_number++;
}
close IN;
$sd=($whole/$number)**0.5;
$sd/=$average;

open(OUT,">quality.txt") or die "Can't open the output file: quality.txt!\n";
print OUT "#the quality distributions' probability#\n";
for($i=0;$i<$max;$i++)
{
	print OUT "$i-$ARGV[1]_rand=";
	for($j=0;$j<$step-1;$j++)
	{
		printf OUT ("%0.5f:",$quality{$i*$step+$j});
	}
	printf OUT ("%0.5f\n",$quality{$i*$step+$j});
}

print OUT "#The reads' average length:#\n";
printf OUT ("length=%0.5f\n",$average);
print OUT "#the sd_ratio(sd/average_length):#\n";
printf OUT ("sd_ratio=%0.5f\n",$sd);

###############################second part
if($sd<=0.001)
{
	close OUT;
	exit;
}

adjust_length();
########function statistics--
sub statistics
{
	my ($char,$position)=@_;
	my $pos;
	if($method==2)  #454
	{
		if($ascii{$char}>82)
		{
			$pos=($position+1)*$step-1;
		}
		else
		{
			$pos=$position*$step+$ascii{$char}-33;
		}
	}
    elsif($method==3)  #pacbio
	{
		if($ascii{$char}>82)
		{
			$pos=($position+1)*$step-1;
		}
		else
		{
			$pos=$position*$step+$ascii{$char}-33;
		}
	}
 	elsif($method==1) #illumina
    {
		if($ascii{$char}>82)
		{
			$pos=($position+1)*$step-1;
		}
		else
		{
			$pos=$position*$step+$ascii{$char}-33;
		}
	}

	if(exists $quality{$pos})
	{
		$quality{$pos}++;
	}
	else
	{
		$quality{$pos}=1;
	}
}

#####function adjust--
sub adjust
{
	my $total=0;
	my $phase=0;
	my $j;
	my ($position)=@_;
	for($j=0;$j<$step;$j++)
	{
		$total+=$quality{$position*$step+$j};
	}

	if($total<=10)  #threshold for max length
	{
		return 0;
	}

	for($j=0;$j<$step;$j++)
	{
		$phase+=$quality{$position*$step+$j};
		$quality{$position*$step+$j}=$phase/$total;
	}
	return 1;
}

###########################rand_length_adjust
sub adjust_length
{
	my ($circle1,$circle2,$sign,$min,$total,$phase);
	if($method==2)
	{
		$min=40;
	}
	elsif($method==1)
	{
		$min=30;
	}
                 elsif($method==3)
	{
		$min=150;
	}
	$total=0;
	$phase=0;
	for($circle1=$max;$circle1>=$min;$circle1--)
	{
		if(exists $rand_length{$circle1})
		{
			$total+=$rand_length{$circle1};
			next;
		}
		$sign=0;
		for($circle2=$circle1-1;$circle2>=$min;$circle2--)
		{
			if(exists $rand_length{$circle2})
			{
				$sign=1;
				last;
			}
		}

		if($sign==0)
		{
			$rand_length{$circle1}=0;
			next;
		}
		$rand_length{$circle1}=($rand_length{$circle1+1}-$rand_length{$circle2})/($circle1+1-$circle2)*($circle1-$circle2)+$rand_length{$circle2};
		$total+=$rand_length{$circle1};
	}
	print OUT "length_rand=";
	for($circle1=$min;$circle1<$max;$circle1++)
	{
		$phase+=$rand_length{$circle1};
		$rand_length{$circle1}=$phase/$total;
		printf OUT ("%0.5f:",$rand_length{$circle1});
	}
	print OUT "1.00000\n";
}

