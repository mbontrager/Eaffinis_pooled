#!/bin/bash
# Default parameters
################################################################
#
# Usage: Place in directory with .bam files.
# bash bwa-stampy_markduplicates.sh NNN REF
# Where NNN is the filename of the bam file (without .bam)
# and REF is the path to the reference genome from which the
# bam files were generated
#
# mbontrager@gmail.com
#
################################################################

RMTMPFILES=true
CALCCOV=true
THREADS=1
RNAME=$1
CURDIR=`pwd`
REF=$2

mkdir -p $RNAME
cd $RNAME

# Sort bam file
samtools sort ../${RNAME}.bam ../${RNAME}-s

# Index sorted bam file
samtools index ../${RNAME}-s.bam

# Mark duplicates and sort
java -Djava.io.tmpdir=`pwd`/tmp -jar \
    /home/lee/bioinformatics/picard-tools-1.135/picard.jar MarkDuplicates \
    INPUT=../${RNAME}-s.bam \
    OUTPUT=${RNAME}-smd.bam \
    METRICS_FILE=${RNAME}-smd.metrics \
    AS=TRUE \
    VALIDATION_STRINGENCY=LENIENT \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    REMOVE_DUPLICATES=TRUE
    TMP_DIR=`pwd`/tmp
samtools sort ${RNAME}-smd.bam ${RNAME}-smds
samtools index ${RNAME}-smds.bam

# Determine Genome Coverage and mean coverage per contig
if $CALCCOV; then
    genomeCoverageBed -ibam ${RNAME}-smds.bam > ${RNAME}-smds.coverage
   # generate table with length and coverage stats per contig (From http://github.com/BinPro/CONCOCT)
    python gen_contig_cov_per_bam_table.py --isbedfiles $REF ${RNAME}-smds.coverage > ${RNAME}-smds.coverage.percontig
fi

# Remove temp files
if $RMTMPFILES; then
    rm ../${RNAME}-s* ${RNAME}-smd.bam 
fi

cd $CURDIR
