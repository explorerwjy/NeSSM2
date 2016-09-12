#usage: perl $0 index-file bwa-file
use strict; use warnings;

my $line;
my $max;
my @array;
my $gi;
my @result;
my %hash;
my $prior="";
my $i=0;
my $j;
my $temp;
my %len;
my %copynumber;
my %name;
my $total;
my $flag;
my $length_test;
my @tt;
my $others;
my $circle;
my $percentage;
my $usage="usage:perl $0 index-file bwa-file";
my $length;

###check input
if(@ARGV!=2)
{
	print "ERROR\n";
	print "$usage\n";
	exit;
}

##read the index file
open(FILE,"<$ARGV[0]") or die "ERROR:Can't open the $ARGV[0] file\n$usage\n";
while(<FILE>)
{
	chomp;
	$line=$_;
	@array=split("\t",$line);
	$len{$array[0]}=$array[1];
	$copynumber{$array[0]}=$array[2];
	$name{$array[0]}=$array[3];
}
close FILE;

####read files
open(FILE,"<$ARGV[1]") or die "ERROR:Can't open the $ARGV[1] file!\n$usage\n";
while(<FILE>)
{
	$line=$_;
	@array=split("\t",$line);
	if(@array<10)   ####don't need line
	{
		next;
	}
	if($array[2] eq "*")
	{
		next;
	}

    ($gi)=$array[2]=~/(\w+\d+\.\d+\.\d+)/;
    #print "$gi\n";
	$hash{$gi}=0.001;

	if((length($prior)!=length($array[9]))||(($prior ne $array[9]) &&($prior ne complement($array[9]))))
	{
		$i++;
		$prior=$array[9];
		$result[$i]=$gi;
		#print "$i\n"
	}
	else
	{
		@tt=split("A",$result[$i]);
		$length_test=@tt;
		$flag=0;
		for($j=0;$j<$length_test;$j++)
		{
			if($gi eq $tt[$j])
			{
				$flag++;
				last;
			}
		}
		if($flag>0)
		{
			next;
		}
		$result[$i]=$result[$i]."A".$gi;
	}

	if($line!~/XA\:Z\:/)
	{
		next;
	}

	($others)=$line=~/XA\:Z\:(.+)\;$/;
	@array=split(";",$others);
	for($circle=0;$circle<@array;$circle++)
	{
		($gi)=$array[$circle]=~/(\w+\d+\.\d+\.\d+)/;
		@tt=split("A",$result[$i]);
                $length_test=@tt;
                $flag=0;
                for($j=0;$j<$length_test;$j++)
                {
                	if($gi eq $tt[$j])
        		{
	                	$flag++;
                                last;
                        }
                }
                if($flag>0)
                {
                	next;
                }
                $result[$i]=$result[$i]."A".$gi;
		$hash{$gi}=0.001;
	}
}
close FILE;
$max=$i;

########handle the data

###unique
for($i=1;$i<=$max;$i++)
{
	if($result[$i]=~/A/)
	{
		next;
	}
	$hash{$result[$i]}+=1;
}
print "Total $i unique reads\n";
###common
#for($i=1;$i<=$max;$i++)
#{
#	if($result[$i]!~/A/)
#	{
#		next;
#	}
#
#	@array=split("A",$result[$i]);
#	$length=@array;
#
#	$total=0;
#	for($j=0;$j<$length;$j++)
#	{
#		$total+=($len{$array[$j]}*$hash{$array[$j]});
#		print "check: $len{$array[$j]} \t $hash{$array[$j]}\n";
#		print "total: $total\n";
#	}

#	for($j=0;$j<$length;$j++)
#	{
#		$hash{$array[$j]}+=$len{$array[$j]}*$hash{$array[$j]}/$total;
#		#print "$total\n";
#	}
#}
#print "Total $i common reads\n";

####delete the plasmid
$total=0;
foreach $gi(keys %hash)
{
	#if($copynumber{$gi}==1)
	#{
	#	delete($hash{$gi});
	#	next;
	#}

	if($hash{$gi}<10)  ###you can change it by yourself
	{
		delete($hash{$gi});
                next;
	}
	$total+=$hash{$gi};
}

foreach $gi (keys %hash)
{
	$percentage=$hash{$gi}/$total;
	if($percentage<0.0001)   ###you can change it by yourself
	{
		delete ($hash{$gi});
	}
}

$total=0;
foreach $gi (keys %hash)
{
	$total+=$hash{$gi};
}

open(OUT,">percentage-16s.txt") or die "Can't create the percentage-16s.txt file!\n";
foreach $gi (sort {$hash{$b}<=>$hash{$a}} keys %hash)
{
	$percentage=$hash{$gi}/$total;
	printf OUT ("%s\t%0.5f\n",$name{$gi},$percentage);
}

#####complement the sequence
sub complement
{
        my ($sequence)=@_;
        my $length;
        my @bases;
        my $b;
        my $a="";

        @bases=split("",$sequence);
        $length=@bases;
        for($b=($length-1);$b>=0;$b--)
        {
                if($bases[$b] eq "A")
                {
                        $a=$a."T";
	        }
                elsif($bases[$b] eq "T")
                {
                        $a=$a."A";
                }
                elsif($bases[$b] eq "C")
                {
                        $a=$a."G";
	        }
                elsif($bases[$b] eq "G")
                {
                        $a=$a."C";
                }
                else
                {
                        $a=$a."N";
                }
        }
        return $a;
}

