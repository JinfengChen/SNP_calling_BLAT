#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"ref:s","qry:s","snp:s","flank:s","blat:s","help");


my $help=<<USAGE;
perl $0 --qry --ref --blat --flank 200
Blat the flanking sequence of SNPs to ancestor genome. Get ancestor state of the SNPs in ancestor.
Flanking sequence store 100+1SNP+100=201, or 500+1SNP+500=1001.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{qry} ||= "../input/HEG4_EG4_A119_A123_NB_SNPs.noRepeats.selectedSNPs.flanking_200.fasta";
$opt{ref} ||= "/rhome/cjinfeng/BigData/00.RD/Transposon_Oryza/OGE_genomes/O.rufipogon/Oryza_rufipogon_W1943_scaffolds_v1.0.fa";
$opt{snp} ||= "../input/HEG4_EG4_A119_A123_NB_SNPs.noRepeats.selectedSNPs.flanking_200.allele";
$opt{blat} ||= "../input/HEG4_EG4_A119_A123_NB_SNPs.noRepeats.selectedSNPs.flanking_200.ORU.best.psl";
$opt{flank} ||= 200; ### flanking length of SNPs, should be half length of previous 

my $refseq=getfastaseq($opt{ref});
my $qryseq=getfastaseq($opt{qry});
#my $allele=allele($opt{snp});
blat2state($opt{blat},$refseq,$qryseq,$opt{flank}/2);
#blat2vcf($opt{blat},$refseq,$qryseq,$allele);

############
#201     0       0       0       0       0       0       0       +       TBGI000006      201     0       201     Chr1    43270923        14046   14247
#1       201,    0,      14046,
#vcf: Chr1    21547   .       A       T       16728.35
sub blat2vcf
{
my ($file,$ref,$qry,$allele)=@_;
my %hash;
my $pre= $1 if ($file=~/(.*)\.blat$/);
open OUT, ">temp.vcf" or die "$!";
open IN, "$file" or die "$!";
for(my $i=0;$i<5;$i++){
   <IN>; ### skip 5 head lines
}
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    next if (exists $hash{$unit[9]}); ### only count for the best hit
    if ($unit[0] == $unit[10]){ ###block size equal qry size, perfect match including SNP
       $hash{$unit[9]}=1;
       my $refallele=substr($ref->{$unit[13]},$unit[15]+$opt{len},1);
       my $qryallele=substr($qry->{$unit[9]},$unit[11]+$opt{len},1);
       print "$unit[13]\t$unit[15]\t.\t$refallele\t$qryallele\t$unit[9]\t$allele->{$unit[9]}->[0]\t$allele->{$unit[9]}->[1]\n";
       my $pos=$unit[15]+$opt{len}+1;
       ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  MSU7_BWA/HEG4_2.3.MSU7_BWA.Chr1.bam 
       print OUT "$unit[13]\t$pos\t.\t$allele->{$unit[9]}->[0]\t$allele->{$unit[9]}->[1]\t100\t.\tAC=8;AF=1.00;AN=8;DP=463\tGT:AD:DP:GQ:PL\t1/1:0,130:130:99:4480,325,0\n";
    }
}
close IN;
close OUT;
`sort -k1,1 -k2,2n temp.vcf > temp.sort.vcf`;
`cat head.vcf temp.sort.vcf > $pre.vcf`;
return \%hash;
}

#S01000000203T   T       C
sub allele
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=[$unit[1],$unit[2]];
}
close IN;
return \%hash;
}


##TBGI000006      Chr01   10949   T       T       T       
sub allele1
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my %al;
    for(my $i=4;$i<@unit;$i++){
       next if ($unit[$i] =~/n/i or $unit[$i] eq $unit[3]);
       $al{$unit[$i]}=1;
    }
    my @snp;
    if (keys %al == 1){
       foreach(keys %al){
           push @snp,$_;
       }
    }
    $hash{$unit[0]}=[$unit[3],$snp[0]];
}
close IN;
return \%hash;
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

sub revcom
{
my ($seq)=@_;
my $rev=reverse $seq;
$rev=~tr/ATGCatgc/TACGtacg/;
return $rev;
}


##read blat file, store ancestor state in hash: SNPs->Species->Allele
##198     1       0       0       0       0       0       0       -       S01000031071A   201     2       201     scaffold49|size484596   484596  442629  442828  1       199,    0,      442629,
##181     0       0       0       0       0       0       0       +       S01000466488T   201     20      201     scaffold27633|size286   286     0       181     1       181,    20,     0,
sub blat2state
{
my ($file,$ref,$qry,$flank)=@_;
my %hash;
my $pre= $1 if ($file=~/(.*)\.best.psl$/);
#`cp head.vcf $pre.genotype.vcf`;
open OUT, ">$pre.temp.state" or die "$!";
open OUT1, ">>$pre.genotype.vcf" or die "$!";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my $line = $_;
    $line =~ s/\,//g;
    my @unit=split("\t",$line);        
    if ($unit[0]/$unit[10] > 0.9 and $unit[17] == 1){ ###block size larger than 90% of qry size and single block match
       my ($chr, $position);
       if ($unit[9]=~/S(\d{2})(\d{9})/){
          $chr = 'scaffold_'.int($1);
          $position = int($2);
       }
       #my $qrymatch1  = $qry->{$unit[9]};
       #my $qrymatch = $unit[8] eq "+" ? substr($qry->{$unit[9]}, $unit[19], $unit{18}) : substr($qry->{$unit[9]}, $unit[19], $unit{18});
       my $qrymatch = $unit[8] eq "+" ? substr($qry->{$unit[9]}, $unit[19], $unit[18]) : substr(revcom($qry->{$unit[9]}), $unit[19], $unit[18]);
       my $refmatch = substr($ref->{$unit[13]}, $unit[20], $unit[18]);
       my $pos = $flank-$unit[19]; ### SNP position in match fragment
       my $mark= 'X' x $pos; 
       my $qryallele = 'N';
       my $refallele = 'N';
       if ($unit[8] eq "+"){
           $qryallele = substr($qrymatch, $pos, 1);
           $refallele = substr($refmatch, $pos, 1);
       }else{
           $qryallele = revcom(substr($qrymatch, $pos, 1));
           $refallele = revcom(substr($refmatch, $pos, 1));
       }
       print OUT ">$_\n";
       print OUT "MAK: $mark\n";
       print OUT "QRY: $qrymatch\n$qryallele\n";
       print OUT "REF: $refmatch\n$refallele\n";

       #print "$unit[13]\t$unit[15]\t.\t$refallele\t$qryallele\t$unit[9]\t$allele->{$unit[9]}->[0]\t$allele->{$unit[9]}->[1]\n";
       #my $pos=$unit[15]+$opt{len}+1;
       ##REF allele is reference allele for SNPs, ALT allele is ancestor allele from the above step
       ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  MSU7_BWA/HEG4_2.3.MSU7_BWA.Chr1.bam
       print OUT1 "$chr\t$position\t.\t$qryallele\t$refallele\t100\t.\tAC=8;AF=1.00;AN=8;DP=463\tGT:AD:DP:GQ:PL\t1/1:0,130:130:99:4480,325,0\n"; 
    }
}
close IN;
close OUT;
close OUT1;
}

