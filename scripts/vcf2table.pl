#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"vcf:s","help");


my $help=<<USAGE;
perl $0 --vcf
Convert vcf to union formot table, so we can write a script to deal with.
INPUT:
Chr1    31071   .       A       G       9183.67 PASS    AC=8;
OUTPUT:
Chr1	31071	S0100031071A	A	G	S0100000092A
USAGE


if ($opt{help} or !$opt{vcf}){
    print "$help\n";
    exit();
}


readtable($opt{vcf});

#Chr1    31071   .       A       G       9183.67 PASS    AC=8;
#reference chromosome: chromosome01
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    my $line=$_;
    next unless ($line=~/^scaffold/);
    my @unit=split("\t",$_);
    my $chr=$unit[0];
    my $c;
    if ($unit[0]=~/(\d+)/){
       $c=sprintf("%02d",$1);
    }
    my $pos=$unit[1];
    my $id =sprintf("%09d",$pos);
    my $ref = $unit[3];
    my $qry = $unit[4];
    $id ="S".$c.$id;
    print "$chr\t$pos\t$id\t$ref\t$qry\t$id\n";
}
close IN;
return \%hash;
}


