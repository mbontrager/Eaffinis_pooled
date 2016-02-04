#!/usr/bin/python

import os, sys, getopt, glob, subprocess, shutil
from Bio import SeqIO

############################################################
# Align de novo assemblies against the reference copepod genome
# and process the output into a format for R analysis
#
# Usage: Run from script directory:
# python -d blast_database -p path_to_assemblies -contig_length
#
# Author: Martin Bontrager
############################################################

def main():
    path = ''
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],'hp:d:l:',['path=',
                                            'database=', 'length='])
    except getopt.GetoptError:
        print('Fix your input')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('-d blast_db -p path_to_assemblies -l min_contig_length')
            sys.exit()
        elif opt in ("-p", "--path"):
            path = str(arg)
        elif opt in ("-d", "--database"):
            db = arg
        elif opt in ("-l", "--length"):
            l = int(arg)
            
    os.chdir(path)
    run('mkdir blast_output')
    sample_list = glob.glob('*.fna')
    
    for i in sample_list:
        filter_by_length(i, l)
        blast(i, db)

# Filter assemblies to only those contigs greater than n base pairs
def filter_by_length(input_file, length):
    infile = open(input_file, 'rU')
    ofname = (input_file.replace('__FINAL_ASSEMBLY.consolidatedContigs.fna',
                           '_3000.fna'))
    out = open(ofname, 'w')

    for i in SeqIO.parse(infile, 'fasta'):
        sequence = i.seq
        if len(sequence) > length:
            SeqIO.write(i, out, 'fasta')

    infile.close()
    out.close()

def blast(input_file, blast_db):
    iname = input_file.replace('__FINAL_ASSEMBLY.consolidatedContigs.fna',
                                '_3000.fna')
    ofname= 'blast_output/' + iname.replace('.fna', '.b6')
    cmd = ('blastn -query ' + iname + ' -db ' + blast_db + ' -out ' + ofname +
           ' -outfmt 6 -max_target_seqs 1 -reward 3 -penalty -4 -gapopen 6 ' +
           '-gapextend 2 -evalue 0.000001')
    run(cmd)
        
# Simplify running bash commands
def run(cmd):
    print('Running command: \n' + cmd + '\n\n') 
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()
    
if __name__ == "__main__":
    main()

    
