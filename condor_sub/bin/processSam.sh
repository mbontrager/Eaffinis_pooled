#!/bin/bash
samtools/samtools view -f 4 /mnt/gluster/bontrager2/bwa/SCE.sam > /mnt/gluster/bontrager2/bwa/SCE_meta.sam
grep -v ^@ /mnt/gluster/bontrager2/bwa/SCE_meta.sam | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > /mnt/gluster/bontrager2/bwa/SCE_meta_1.fastq
grep -v ^@ /mnt/gluster/bontrager2/bwa/SCE_meta.sam | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > /mnt/gluster/bontrager2/bwa/SCE_meta_2.fastq
samtools/samtools view -h -F 4 /mnt/gluster/bontrager2/bwa/SCE.sam > /mnt/gluster/bontrager2/bwa/SCE_Eaff.sam
samtools/samtools view -bS /mnt/gluster/bontrager2/bwa/SCE_Eaff.sam -o /mnt/gluster/bontrager2/bwa/SCE_Eaff.bam
samtools/samtools sort /mnt/gluster/bontrager2/bwa/SCE_Eaff.bam /mnt/gluster/bontrager2/bwa/SCE_Eaff.sorted
