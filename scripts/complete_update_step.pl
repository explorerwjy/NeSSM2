#!/usr/bin/perl -w
use strict;

my $usage =  "usage: perl complete_updata_step.pl whole_pathway_datebase_dir/";

if(@ARGV != 1) {print $usage."\n"; exit(0);}

my $db_dir=$ARGV[0];

my $temp=(); my $genome=();

if($db_dir=~/\/$/) {$temp=$db_dir."temp/";$genome=$db_dir."genome/";}
else {$db_dir.="/"; $temp=$db_dir."temp/";$genome=$db_dir."genome/";}

print $db_dir."\n";

if(-e $temp) {system("rm -rf $temp");}
system("mkdir $temp");

system("wget -c -t 0 ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.fna.tar.gz -P $db_dir");

my $old=(); my $new_1=(); my $new_2=();
$old=$db_dir."all.fna.tar.gz";  #print $old."\n";
$new_1=$db_dir."bacteria.fna.tar.gz";  #print $new_1."\n";

system("tar -zxvf $old -C $temp");

system("mv $old $new_1");

system("wget -c -t 0 ftp://ftp.ncbi.nih.gov/genomes/Viruses/all.fna.tar.gz -P $db_dir");

$new_2=$db_dir."viruses.fna.tar.gz";

system("tar -zxvf $old -C $temp");

system("mv $old $new_2");


my $cmd=q{ftp -n<<!
open ftp.ncbi.nih.gov
user anonymous nuaaecust@163.com
prompt
ls /genomes/Fungi/ list_fungi
bye
!};

system($cmd);

system("touch ./signal");

my $file="signal";

if(-f $file) {system("perl dl_fungi.pl ./list_fungi $temp");}
else {system("echo 'No signal file, please manual download Fungi sequence data.'");}

system("rm $file");


if(-e $genome) {system("rm -rf $genome");}
system("mv $temp $genome");

system("perl mk_index.pl $genome");

system("rm $new_1");
system("rm $new_2");

system("echo 'Last update: ' > log");
$cmd=q{echo $(date +%Y"/"%m"/"%d" "%k":"%M":"%S) >> log};
system($cmd);


