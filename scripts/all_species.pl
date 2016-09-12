#!/usr/bin/perl -w
use strict;

my $usage =  "usage: perl mk_index.pl whole_pathway_datebase_dir/ output_dir/";
my $index;
my $index_new;
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
	$index=$ARGV[1]."all_species";
    	$index_new=$ARGV[1]."all_species.new";
}
else
{
    	$index=$ARGV[1]."/"."all_species";
    	$index_new=$ARGV[1]."/"."all_species.new";
}

open(INDEX,">$index_new");

if($ARGV[0]!~/\/$/)
{
    	$path=$ARGV[0]."/";
}
else
{
    	$path = $ARGV[0];
}

scan_dir($path);

close(INDEX);

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

	    		scan_dir($list[$i]);
		}
		elsif(-f $list[$i])
		{
	    		if($list[$i]=~/fna$/)
			{
				print INDEX "$list[$i]\n";
				#read_fna($list[$i],$number);
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
			($gi)=$fna[$x]=~/gi\|(\d+)\|ref/;
	    		my @split_1=split /\,/,$fna[$x];
	    		my @split_2=split /\|\s/,$split_1[0];
	    		$tmp_name=$split_2[1];
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
    	print INDEX "$gi\t$length\t$_[1]\t$tmp_name\t$_[0]\n";

    	close FNA;
}
