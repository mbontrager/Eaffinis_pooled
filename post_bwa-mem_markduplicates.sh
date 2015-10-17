#!/bin/bash
# Default parameters
RMTMPFILES=true
CALCCOV=true
THREADS=1
RNAME=$1
QNAME="Eaff"
CURDIR=`pwd`

mkdir -p $RNAME
cd $RNAME

# Index sorted bam file
samtools index ../${RNAME}_${QNAME}.sorted.bam

# Mark duplicates and sort
java -Djava.io.tmpdir=`pwd`/tmp -jar \
    /home/lee/bioinformatics/picard-tools-1.135/picard.jar MarkDuplicates \
    INPUT=../${RNAME}_${QNAME}.sorted.bam \
    OUTPUT=${RNAME}_${QNAME}-smd.bam \
    METRICS_FILE=${RNAME}_${QNAME}-smd.metrics \
    AS=TRUE \
    VALIDATION_STRINGENCY=LENIENT \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    REMOVE_DUPLICATES=TRUE
    TMP_DIR=`pwd`/tmp
samtools sort ${RNAME}_${QNAME}-smd.bam ${RNAME}_${QNAME}-smds
samtools index ${RNAME}_${QNAME}-smds.bam

# Determine Genome Coverage and mean coverage per contig
if $CALCCOV; then
    genomeCoverageBed -ibam ${RNAME}_${QNAME}-smds.bam > ${RNAME}_${QNAME}-smds.coverage
    awk 'BEGIN {pc=""}
    {
        c=$1;
        if (c == pc) {
            cov=cov+$2*$5;
        } else {
            print pc,cov;
            cov=$2*$5;
        pc=c}
    } END {print pc,cov}' ${RNAME}_${QNAME}-smds.coverage | tail -n +2 > ${RNAME}_${QNAME}-smds.coverage.percontig
fi

# Remove temp files
if $RMTMPFILES; then
    rm ../${RNAME}_${QNAME}.sorted.bam \
       ${RNAME}_${QNAME}-smd.bam \
fi

cd $CURDIR
