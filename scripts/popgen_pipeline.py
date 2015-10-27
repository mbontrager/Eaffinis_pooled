#!/usr/bin/python

import sys, getopt, glob, subprocess, shutil
from optparse import OptionParser

################################################################################
# Wrap script for pooled population genomics analysis
#
# Author: Martin Bontrager
#
# Usage: python popgen_pipeline.py -r /full/path/to/reference/genome.fa -d
#                                  /full/path/to/bamfile/dir/
#
################################################################################

def main():
    
    parser = OptionParser(usage='python popgen_pipeline.py -r' +
                          ' <path_to_ref_genome> -d <path_to_bam_directory>')
    parser.add_option('-r', '--reference', 
                        dest='ref',
                        help='foo help')
    parser.add_option('-d', '--bamdir',
                      dest='dir',
                      help='foo help')
    (options, args) = parser.parse_args()

    if not options.ref:   # if ref is not given
        parser.error('Reference genome not given')
        sys.exit(2)
    if not options.dir:   # if dir is not given
        parser.error('Bam directory not give')
        sys.exit(2)

    options.dir.rstrip('/')
    output = glob.glob(options.dir + "/*.bam")

    sample_list = []
    for f in output:
        sample_list.append(f.replace(options.dir, '').replace('.bam', ''))

    #mark_duplicates(sample_list, options.ref, options.dir)
    #filter_unmapped(sample_list, options.dir)
    count_unmapped(sample_list, options.dir)

# Simplify running bash commands
def run(cmd):
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()

def mark_duplicates(sample_list, ref, dir):
    for n in sample_list:
        cmd = ('bash bwa-stampy_markduplicates.sh ' + dir + ' ' + n )
        run(cmd)
        cmd = ('python gen_contig_cov_per_bam_table.py --isbedfiles ' +
               ref + ' ' + dir + n + '/' + n + '-smds.coverage > ' + dir + n +
               '/' + n + '-smds.coverage.py.percontig')
        run(cmd)

def filter_unmapped(sample_list, dir):
    for i in sample_list:

        cmd = ('samtools view -b -f 4 ' + dir + i + '.bam -o ' + dir +
               'unmapped/' + i + '_unmapped.bam')
        run(cmd)
        cmd = ('samtools sort ' + dir + 'unmapped/' + i + '_unmapped.bam ' +
               dir + 'unmapped/' + i + '-s')
        run(cmd)
        cmd = ('samtools index ' + dir + 'unmapped/' + i + '-s.bam')
        run(cmd)
        cmd = ('rm ' + dir + 'unmapped/' + i + '_unmapped.bam')
        run(cmd)

def count_unmapped(sample_list, dir):
    for i in sample_list:
        cmd = ('samtools idxstats ' + dir + i + '/' + i + '-smds.bam > ' + dir +
               i + '/' + i + '-smds.mappingstats.tsv')
        run(cmd)

    
if __name__ == "__main__":
    main()
