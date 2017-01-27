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
python SNP_call_by_BLAT.py --ref Cclementina_v1.0_scaffolds.fa --ass citrus.contigs.fasta --vcf dbSNP.clean.Het.vcf --flank 200 --output ass.genome.SNP

We genotype dbSNP in a new assembly by mapping the flanking sequence of SNP to new assembly using BLAT.
1. convert vcf into table
2. extract flanking sequence from reference
3. BLAT map flanking sequence to new genome
4. extact genotype from new genome and write to vcf
    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 60 --lines %s --interval 120 --task 1 --mem 15G --time 100:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

def write_file(cmd, filename):
    ofile = open(filename, 'w')
    print >> ofile, cmd
    ofile.close()

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', '--ref')
    parser.add_argument('--assembly', '--asm')
    parser.add_argument('--vcf')
    parser.add_argument('--flank', help='set 200 bp flanking will extract 100 bp from each side')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.reference) > 0 and len(args.assembly) > 0 and len(args.vcf) > 0
    except:
        usage()
        sys.exit(2)
 
    if not args.output:
        args.output = 'BLAT.SNP'

    if not args.flank:
        args.flank  = 200


    prefix = args.output
    args.output    = os.path.abspath(args.output)
    args.vcf       = os.path.abspath(args.vcf)
    args.reference = os.path.abspath(args.reference)
    args.assembly  = os.path.abspath(args.assembly)

    fatherpath=os.path.dirname(os.path.abspath(sys.argv[0]))
    script='%s/scripts' %(fatherpath)

    #convert vcf into table
    step1_shell = '%s.step1.sh' %(args.output)
    cmd = 'perl %s/vcf2table.pl --vcf %s > %s.table' %(script, args.vcf, args.output)
    write_file(cmd, step1_shell)
    if not os.path.exists('%s.table' %(args.output)):
        runjob(step1_shell, 1)
 
    #extract flanking sequence of SNP
    step2_shell = '%s.step2.sh' %(args.output)
    cmd = 'perl %s/SNP_flanking.pl --table %s.table --ref %s --flank %s' %(script, args.output, args.reference, args.flank)
    write_file(cmd, step2_shell)
    if not os.path.exists('%s.flanking_%s.fasta' %(args.output, args.flank)):
        runjob(step2_shell, 1)
 
    #blat flanking sequence to new assembly
    step3_shell = '%s.step3.sh' %(args.output)
    os.system('perl ~/BigData/software/bin/fastaDeal.pl --cuts 10000 %s.flanking_%s.fasta' %(args.output, args.flank))
    subfiles = glob.glob('%s.flanking_%s.fasta.cut/*.fasta.*' %(args.output, args.flank))
    cmds = []
    for fa in sorted(subfiles):
        cmds.append('/opt/linux/centos/7.x/x86_64/pkgs/blat/35/bin/blat -noHead -minIdentity=95 %s %s %s.psl' %(args.assembly, fa, fa))
        cmds.append('perl %s/bestAlign.pl --cutoff 0.90 %s.psl > %s.best.psl' %(script, fa, fa))
    cmd = '\n'.join(cmds)
    write_file(cmd, step3_shell)
    if not os.path.exists('%s.best.psl' %(args.output)):
        runjob(step3_shell, 2)
        os.system('cat %s.flanking_%s.fasta.cut/*.best.psl > %s.best.psl' %(args.output, args.flank, args.output))


    #extract SNP and write vcf
    step4_shell = '%s.step4.sh' %(args.output) 
    cmd = 'perl %s/blat2state.pl --qry %s.flanking_%s.fasta --ref %s --blat %s.best.psl --flank %s' %(script, args.output, args.flank, args.assembly, args.output, args.flank)
    write_file(cmd, step4_shell)
    if not os.path.exists('%s.genotype.vcf' %(args.output)):
        runjob(step4_shell, 1)

if __name__ == '__main__':
    main()

