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

    # Add all .bam files in the input directory to the samples list
    sample_list = []
    for f in output:
        sample_list.append(f.replace(options.dir, '').replace('.bam', ''))
    
    # Modify this path to point at your picard .jar file 
    picard_path = '/home/lee/bioinformatics/picard-tools-1.135/picard.jar'
    
    # Mark duplicates in the bam files, and sort, index the files   
    #mark_duplicates(sample_list, options.ref, options.dir, picard_path)
    
    # Filter unmapped reads from the bam file. And perform counts of unmapped
    # reads to understand mapping patterns
    #filter_unmapped(sample_list, options.dir)
    #count_unmapped(sample_list, options.dir)
    
    # Add read group headers to bam files. This is a very naive implementation
    # of adding read groups. I don't use the Illumina run information. As such,
    # there are aspects of read correction in GATK that may not work well.
    #add_headers(sample_list, options.dir)
    
    # Find indel targets for realignment and then realign around them.
    #target_creator(sample_list, options.ref, options.dir)
    #realign(sample_list, options.ref, options.dir)
    
    # Recalibrate base quality scores
    #recalibrate_qual(sample_list, options.ref, options.dir)
    
    #Call variants
    call_snps(sample_list, options.ref, options.dir)

# Simplify running bash commands
def run(cmd):
    print('Running command: \n' + cmd + '\n\n') 
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()

# Samtools sort bam file
def sort(inDir, outDir, f, h='-s'):
    s = f.replace('.bam', '')
    cmd = ('samtools sort ' + inDir + f + ' ' + outDir + s + h)
    run(cmd)

# Samtools index bam file
def index(filepath):
    cmd = ('samtools index ' + filepath)
    run(cmd)    

#Sort, index, run Picard MarkDuplicats, then sort and index again
# Also calculate coverage per contig    
def mark_duplicates(sample_list, ref, dir, picard_path):
    for n in sample_list:
        run('mkdir ' + dir + n)
        sort(dir, (dir + n + '/'), (n + '.bam'))
        index(dir + n + '/' + n + '-s.bam')
        
        cmd = ('java -Djava.io.tmpdir=' + dir + n + '/tmp -jar ' + 
        picard_path + ' MarkDuplicates ' +
        'INPUT=' + dir + n + '/' + n + '-s.bam ' +
        'OUTPUT=' + dir + n + '/' + n + '-smd.bam ' +
        'METRICS_FILE=' + dir + n + '/' + n + '-smd.metrics ' +
        'AS=TRUE ' +
        'VALIDATION_STRINGENCY=LENIENT ' + 
        'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 ' + 
        'REMOVE_DUPLICATES=TRUE ' + 
        'TMP_DIR=' + dir + n + '/tmp')
        run(cmd)
      
        odir = (dir + n + '/')
        smd = n + '-smd.bam'
        smds = odir + n + '-smds.bam'
        
        sort(odir, odir, smd, h='s')
        index(smds)
        
        # Calculate coverage
        cmd = ('genomeCoverageBed -ibam ' + smds + ' > ' + odir + n +
               '-smds.coverage')
        run(cmd)

        # Remove intermediate files:
        run('rm ' + odir + smd)
        run('rm ' + odir + n + '-s.bam.bai')
        run('rm ' + odir + n + '-s.bam')

# Use samtools to filter unmapped reads from a bam file        
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

# Gather stats on unmapped reads        
def count_unmapped(sample_list, dir):
    for i in sample_list:
        cmd = ('samtools idxstats ' + dir + i + '/' + i + '-smds.bam > ' + 
                dir + i + '/' + i + '-smds.mappingstats.tsv')
        run(cmd)

# Find indels in an alignment
def target_creator(sample_list, ref, dir):
    for i in sample_list:
        
        cmd = ('samtools index ' + dir + i + '/' + i + '-smds_withRG.bam')
        run(cmd)
        
        cmd = ('java -jar $GATK ' + 
               '-T RealignerTargetCreator ' +
               '-R ' + ref + ' ' +
               '-I ' + dir + i + '/' + i + '-smds_withRG.bam ' +
               '-o ' + dir + i + '/' + i + '_realigned.list')
        run(cmd)

# Add read group headers to a bam file **NAIVE IMPLEMENTATION**        
def add_headers(sample_list, dir):
    for i in sample_list:
        file = dir + i + '/' + i + '-smds.bam'
        cmd = ('java -jar  ' + 
        '/home/lee/bioinformatics/picard-tools-1.135/picard.jar ' + 
        'AddOrReplaceReadGroups ' +
        'I=' + file + ' ' +
        'O=' + dir + i + '/' + i + '-smds_withRG.bam ' +
        'SORT_ORDER=coordinate ' + 
        'RGID=' + i + ' ' + 
        'RGLB=' + i + '1 ' +  
        'RGPL=ILLUMINA ' + 
        'RGPU=' + i + 'seq ' + 
        'RGSM=' + i)
        run(cmd)
        cmd = ('rm ' + file)
        run(cmd)

# Realign reads around indel targets
def realign(sample_list, ref, dir):
    for i in sample_list:
        file = dir + i + '/' + i + '-smds_withRG.bam'
        relist = dir + i + '/' + i + '_realigned.list'
        cmd = ('java -jar $GATK ' + 
               '-T IndelRealigner ' +
               '-R ' + ref + ' ' +
               '-I ' + file + ' ' +
               '-targetIntervals ' + relist + ' ' +
               '-o ' + dir + i + '/' + i + '_realigned.bam')
        run(cmd)
        cmd = ('rm ' + file)
        run(cmd)
 
# Recalibrate base quality scores. Unused since we don't know polymorphic sites
def recalibrate_qual(sample_list, ref, dir):
    for i in sample_list:
        file = dir + i + '/' + i +'_realigned.bam'
        cmd = ('java -jar $GATK ' + 
        '-T BaseRecalibrator ' +
        '-R ' + ref + ' ' + 
        '-I ' + file + ' ' + 
        '-o ' + dir + i + '/' + i + '_recal_data.table')
        
        print(cmd)
        
        cmd = ('java -jar $GATK ' +
        '-T PrintReads ' + 
        '-R ' + ref + ' ' +
        '-I ' + file + ' ' + 
        '-BQSR ' + dir + i + '/' + i + '_recal_data.table ' + 
        '-o ' + dir + i + '/' + i + '_recal.bam')
        print(cmd)
        
#First round of SNP calling with HaplotypeCaller
def call_snps(sample_list, ref, dir):
    for i in sample_list:
        file = dir + i + '/' + i +'_realigned.bam'
        cmd = ('java -jar $GATK ' + 
        '-T HaplotypeCaller ' + 
        '-R ' + ref + ' ' + 
        '-I ' + file + ' ' + 
        '-ploidy 200 ' + 
        '--genotyping_mode DISCOVERY ' + 
        '-stand_emit_conf 10 '  + 
        '-stand_call_conf 30 ' + 
        '-o ' + dir + i + '/' + i + '_variants.vcf')
        print(cmd)

if __name__ == "__main__":
    main()
