#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"table:s","ref:s","help");


my $help=<<USAGE;
perl $0 --table --ref --flank
Read snp table and ref genome, generate a fasta file with 100+1SNP+100 sequence and 500+1SNP+500.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{table} ||="../input/HEG4_EG4_A119_A123_NB_SNPs.noRepeats.selectedSNPs.table";
$opt{ref} ||="/rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa";
$opt{flank} ||= '200';
my $refseq=getfastaseq($opt{ref});
readtable($opt{table},$refseq,$opt{flank});

#reference chromosome: chromosome01
#chromosome01    92      S01000000092    A       T
#position 92 is 1-based, so 92-1=91, 0-based, is used in substr as start position.
sub readtable
{
my ($file,$refseq, $flanklen)=@_;
my $pre=$1 if ($file=~/(.*)\.table/);
open OUT1, ">$pre.flanking_$flanklen.fasta" or die "$!";
open OUT2, ">$pre.flanking_$flanklen.allele" or die "$!";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $base  =substr($refseq->{$unit[0]},$unit[1]-1,1);
    next if ($unit[1]- ($flanklen/2 + 1) < 0); ### not have 100 before SNPs, skip these cases
    my $flank1=substr($refseq->{$unit[0]},$unit[1]- ($flanklen/2 + 1), $flanklen+1);
    $flank1=formatseq($flank1,100);
    next if (length $flank1 < ($flanklen + 1)); ### not have 100 after SNPs, skip these cases
    #print "$_\t$base\n"; 
    print OUT1 ">".$unit[2].$base."\n$flank1\n";
    
    ##allele
    my $al= $base eq $unit[3] ? $unit[4] : $unit[3];
    print OUT2 $unit[2].$base."\t$base\t$al\n";
}
close IN;
close OUT1;
close OUT2;
}
 
sub getfastaseq
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}=$seq;
}
$/="\n";
return \%hash;
}


sub formatseq
{
### format a single line sequence into lines with user specific length
my ($seq,$step)=@_;
my $length=length $seq;
my $run=int ($length/$step);
my $newseq;
for(my $i=0;$i<=$run;$i++){
   my $start=$i*$step;
   my $line=substr($seq,$start,$step);
   $newseq.="$line\n";
}
return $newseq;
}


