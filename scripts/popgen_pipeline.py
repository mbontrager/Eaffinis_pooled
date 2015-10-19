#!/usr/bin/python

import os, sys, getopt, glob, subprocess, shutil

############################################################
# Wrap script for pooled population genomics analysis
#
# Author: Martin Bontrager
#
############################################################

def main():
    path = os.path.dirname(os.path.realpath(__file__))
    ref = '/media/lee/new_york/EAFFGENOME/assembl/Eaff_11172013.genome.fa'

    sample_list = [name for name in os.listdir(path)
            if os.path.isdir(os.path.join(path, name))]

# Simplify running bash commands
def run(cmd):
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()

def mark_duplicates(sample_list):
    for n in sample_list:
        #cmd = ('bash bwa-stampy_markduplicates.sh ' + n)
        #run(cmd)
        cmd = ('python gen_contig_cov_per_bam_table.py --isbedfiles ' +
               ref + ' ' + n + '/' + n + '-smds.coverage > ' + n + '/' +
               n + '-smds.coverage.py.percontig')

    
if __name__ == "__main__":
    main()
