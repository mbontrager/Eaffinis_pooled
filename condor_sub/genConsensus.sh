#!/bin/bash
samtools/samtools mpileup -uf /mnt/gluster/bontrager2/eaff_contigs/Eaff_11172013.genome.fa /mnt/gluster/bontrager2/CBE_realigned.bam | bcftools view -cg - | vcfutils.pl vcf2fq > /mnt/gluster/bontrager2/CBE_cns.fq
