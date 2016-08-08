#!/bin/bash
bwa/bwa mem -M -t 20  /mnt/gluster/bontrager2/eaff_contigs/Eaff_11172013.genome.fa /mnt/gluster/bontrager2/bwa/MIE_R1.fastq.gz /mnt/gluster/bontrager2/bwa/MIE_R2.fastq.gz > MIE.sam

