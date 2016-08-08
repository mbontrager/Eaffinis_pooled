#!/bin/bash
bwa/bwa mem -M -t 32 \
/mnt/gluster/bontrager2/eaff_masked/Eaff_11172013.genome.masked.fa \
/mnt/gluster/bontrager2/popgen_samples/$1/$1_R1.trimmed.paired.fq.gz \
/mnt/gluster/bontrager2/popgen_samples/$1/$1_R2.trimmed.paired.fq.gz \
> /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.sam

#Process file
samtools/samtools view -bS /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.sam \
-o /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.bam
rm /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.sam
rm /mnt/gluster/bontrager2/popgen_samples/$1/$1_R*