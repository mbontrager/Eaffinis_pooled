#!/bin/bash
samtools/samtools view -bS /mnt/gluster/bontrager2/bwa/MAE_Eaff.sam -o /mnt/gluster/bontrager2/bwa/MAE_Eaff.bam
samtools/samtools sort /mnt/gluster/bontrager2/bwa/MAE_Eaff.bam /mnt/gluster/bontrager2/bwa/MAE_Eaff.sorted
