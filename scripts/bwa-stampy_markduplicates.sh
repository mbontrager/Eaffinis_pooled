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
RNAME=$2
DNAME=$1

mkdir -p ${DNAME}${RNAME}

# Sort bam file
samtools sort ${DNAME}${RNAME}.bam ${DNAME}${RNAME}/${RNAME}-s

# Index sorted bam file
samtools index ${DNAME}${RNAME}/${RNAME}-s.bam

# Mark duplicates and sort
java -Djava.io.tmpdir=`pwd`/tmp -jar \
    /home/lee/bioinformatics/picard-tools-1.135/picard.jar MarkDuplicates \
    INPUT=${DNAME}${RNAME}/${RNAME}-s.bam
    OUTPUT=${DNAME}${RNAME}/${RNAME}-smd.bam \
    METRICS_FILE=${DNAME}${RNAME}/${RNAME}-smd.metrics \
    AS=TRUE \
    VALIDATION_STRINGENCY=LENIENT \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    REMOVE_DUPLICATES=TRUE
    TMP_DIR=`pwd`/tmp
samtools sort ${DNAME}${RNAME}/${RNAME}-smd.bam ${DNAME}${RNAME}/${RNAME}-smds
samtools index ${DNAME}${RNAME}/${RNAME}-smds.bam

# Determine Genome Coverage and mean coverage per contig
if $CALCCOV; then
    genomeCoverageBed -ibam ${DNAME}${RNAME}/${RNAME}-smds.bam > \
    ${DNAME}${RNAME}/${RNAME}-smds.coverage
fi

# Remove temp files
if $RMTMPFILES; then
    rm ${DNAME}${RNAME}/${RNAME}-s.bam \
       ${DNAME}${RNAME}/${RNAME}-smd.bam \
       ${DNAME}${RNAME}/${RNAME}-s.bai
fi

# Use GATK to realign around indels

