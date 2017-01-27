#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir

def usage():
    test="name"
    message='''
python Summary_Haplotype_SNP_by_Het_dbSNP.py --dbSNP dbSNP.clean.Het.vcf --phasedSNP scaffold_all.snp.list --hapSNP ass.genome.SNP.genotype.vcf

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 60 --lines 2 --interval 120 --task 1 --mem 15G --time 100:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

#SNP	donor	hap1	hap2
#scaffold_1:525279	C/A	C	A
def read_phasedSNP(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith(r'scaffold'): 
                unit = re.split(r'\t',line)
                snp  = unit[0]
                data[snp] = [unit[2], unit[3]]
    return data

#scaffold_1	74012	.	A	T	2903.67	PASS	.	GT:AD:DP:GQ:PL	0/1:4,12:16:99:397,0,108
def readvcf(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith(r'scaffold'): 
                unit = re.split(r'\t',line)
                snp  = '%s:%s' %(unit[0], unit[1])
                data[snp] = [unit[3], unit[4]]
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dbSNP')
    parser.add_argument('--hapSNP')
    parser.add_argument('--phasedSNP')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.dbSNP) > 0
    except:
        usage()
        sys.exit(2)

    dbSNP = readvcf(args.dbSNP)
    hapSNP= readvcf(args.hapSNP)
    phasedSNP = read_phasedSNP(args.phasedSNP)
    dbSNP_n  = len(dbSNP.keys())
    hapSNP_n = len(hapSNP.keys())
    phasedSNP_n = len(phasedSNP.keys())
    compared = 0
    hap1     = 0
    hap2     = 0
    other    = 0
    for site in sorted(dbSNP.keys()):
        if phasedSNP.has_key(site):
            if hapSNP.has_key(site): 
                compared += 1    
                if hapSNP[site][1] == phasedSNP[site][0]:
                    hap1 += 1
                elif hapSNP[site][1] == phasedSNP[site][1]:
                    hap2 += 1
                else:
                    other += 0 
    print 'Compared Phased SNP: %s, %s of all Het SNP and %s of all Phased SNP' %(compared, str(float(compared)/dbSNP_n), str(float(compared)/phasedSNP_n))
    print '%s (%s) is haplotype1 and %s (%s) is haplotype2' %(str(hap1), float(hap1)/compared, str(hap2), str(float(hap2)/compared))
 
if __name__ == '__main__':
    main()

