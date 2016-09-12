#This script is used for add 16s rRNA copy number to the index file for simulation
#usage: perl add_cpn.pl index 16s_rRNA_copy_number_counting_file

#use strict;
use warnings;
#use diagnostics

my $index = $ARGV[0];
my $info = $ARGV[1];
my $index16=$ARGV[0]."cpn";
my $specie = 0;
my $flag = 0;



open($INFILE1, "$ARGV[0]") or die "Cannot open '$INFILE1'";
open(NEW,">$index16");

while($line1 = <$INFILE1>)
{
    chomp($line1);
    @array=split("\t",$line1);
    open($INFILE2, "$ARGV[1]") or die "Cannot open '$INFILE2'";
    while($line2 = <$INFILE2>)
    {
        chomp($line2);
        my @words = split /\t/ , $line2;
        if($line1 =~ /$words[0]/)  #genome name matched
        {
            print NEW "$array[0]\t$array[1]\t$words[1]\t$array[3]\t$array[4]\n";
            $specie++;
            print "$specie specie add in the index file\n";
            $flag = 1;
            close $INFILE2;
            last
        }
    }
    if($flag == 0)
    {
        print NEW "$array[0]\t$array[1]\t1\t$array[3]\t$array[4]\n";
    }
    $flag = 0;

}
print "Total $specie specie add in the index file\n";

close NEW;
close $INFILE1;

sub divide_two
{
    my $subin = shift;
    my @words = split / / , $subin;
    my $first = $words[0];
    splice (@words, 0, 1);

    my $second = "";
    for my $word (@words)
    {
        $second = $second . $word . ' ';
    }
    chomp($second);

    return ($first ,$second);
}

sub detect
{
    my $subin = shift;
    my @words = split /\// , $subin;
    my $tmp = $words[0];
    splice (@words, 0, 1);
    for my $word (@words)
    {
        $tmp = $tmp .'_' . $word;
    }
    chomp ($tmp);
    return ($tmp);
}


