#!/bin/bash
java -Xmx8g -jar /home/lee/bioinformatics/picard-tools-1.135/picard.jar \
AddOrReplaceReadGroups \
I=/media/lee/new_york/population_genomics/2016-04_BWA_stampy_alignments/$1.bwa.stampy.sorted.cleaned.bam \
O=/media/lee/new_york/population_genomics/2016-04_BWA_stampy_alignments/$1.bwa.stampy.sorted.cleaned.RG.bam \
SORT_ORDER=coordinate \
RGID=$1 \
RGLB=$1a \
RGPL=ILLUMINA \
RGPU=$1seq \
RGSM=$1

