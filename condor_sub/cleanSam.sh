#!/bin/bash
samtools/samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 \
-b /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.stampy.sorted.bam \
-o /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.stampy.sorted.cleaned.bam
rm /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.stampy.sorted.bam

