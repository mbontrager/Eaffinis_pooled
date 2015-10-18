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

    sample_list = glob.glob('*.bam')

    for filename in sample_list:
        n = filename.replace('.bam', '')
        cmd = ('bash bwa-stampy_markduplicates.sh ' + n)
        run(cmd)

# Simplify running bash commands
def run(cmd):
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()
    
if __name__ == "__main__":
    main()
