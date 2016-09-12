
#!/usr/bin/perl -w
use strict;

my $usage =  "usage: perl mk_16_index.pl whole_pathway_16s_datebase_dir/ output_dir/";
my $index16;
my $index16_new;
my $line;
my $path="";
my $gi;
my $number=0;
my @array;

if(@ARGV != 2)
{
	print $usage."\n"; exit(0);
}

if($ARGV[1]=~/\/$/)
{
	$index16=$ARGV[1]."index16";
    $index16_new=$ARGV[1]."index16.new";
}
else
{
    $index16=$ARGV[1]."/"."index16";
    $index16_new=$ARGV[1]."/"."index16.new";
}

open(INDEX16,">$index16_new");

if($ARGV[0]!~/\/$/)
{
    	$path=$ARGV[0]."/";
}
else
{
    	$path = $ARGV[0];
}

scan_dir($path);

close(INDEX16);

##adjust
open(NEW,">$index16");
open(OLD,"<$index16_new");
while(<OLD>)
{
	chomp;
	$line=$_;
	@array=split("\t",$line);
	if(($array[3]=~/plasmid/)||($array[3]=~/fragment/))
	{
		print NEW "$array[0]\t$array[1]\t1\t$array[3]\t$array[4]\n";
	}
	else
	{
		print NEW "$array[0]\t$array[1]\t1\t$array[3]\t$array[4]\n";
	}
}
close NEW;
close OLD;
system("rm $index16_new");


sub scan_dir
{
	my $a = $_[0];
    	opendir(DIR,$a);
    	my @dir = readdir(DIR);
    	my @list;
    	my $k=0;
    if($a!~/\/$/)
    {
		$a.="/";
    }

    for(my $i=0;$i<@dir;$i++)
	{
		if($dir[$i]=~/^\./)
		{
	    	next;
		}
		$dir[$i]=$a.$dir[$i];
		$list[$k]=$dir[$i];
		#print "$k $list[$k]\n";
		$k++;
    }

    for(my $i=0;$i<@list;$i++)
	{
		if($list[$i]=~/^\./)
		{
	    	next;
		}
		elsif(-d $list[$i])
		{
	    	#print "oh,another dir!!\n";
	    	scan_dir($list[$i]);
		}
		elsif(-f $list[$i])
		{
	    	if($list[$i]=~/fasta$/)
			{
				#print "$i oh,a fasta file!!\n";
				read_fna($list[$i],$number);
            }
		}
		else
		{
	    		next;
		}
    }

    close DIR;
	$number++;
}


sub read_fna
{
    open(FNA,$_[0]);
    my @fna=<FNA>;
    my $length=0;
    my $tmp_name="";
    for(my $x=0;$x<=$#fna;$x++)
	{
		if($x==0)
		{
	    	chomp $fna[$x];
		    ($gi)=$fna[$x]=~/(w+)(d+)\.(d+)\.(d+)/;
	    	my @split_1=split / /,$fna[$x];
	    	my @split_2=split /\;/,$fna[$x];
	    	$gi = $split_1[0];
  		    $gi=~(/^>([^ ]+)/);
            $gi=$1;
	    	$tmp_name=$split_2[6];
		}
		else
		{
	    	chomp($fna[$x]);
			$length+=length($fna[$x]);
		}
    }
    close(FNA);
    if($tmp_name eq "")
    {
        print STDERR "File $_[0] title is not regular\n";
    }
    print INDEX16 "$gi\t$length\t$_[1]\t$tmp_name\t$_[0]\n";

    close FNA;
}
