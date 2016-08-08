#!/bin/bash        
java -Xmx32g -jar trimmomatic-0.36.jar PE -threads 20 -phred33 \
/mnt/gluster/bontrager2/popgen_samples/$1/$1_R1.fastq.gz \
/mnt/gluster/bontrager2/popgen_samples/$1/$1_R2.fastq.gz \
/mnt/gluster/bontrager2/popgen_samples/$1/$1_R1.trimmed.paired.fq.gz \
/mnt/gluster/bontrager2/popgen_samples/$1/$1_R1.trimmed.unpaired.fq.gz \
/mnt/gluster/bontrager2/popgen_samples/$1/$1_R2.trimmed.paired.fq.gz \
/mnt/gluster/bontrager2/popgen_samples/$1/$1_R2.trimmed.unpaired.fq.gz \
ILLUMINACLIP:/mnt/gluster/bontrager2/popgen_samples/adapters.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

rm /mnt/gluster/bontrager2/popgen_samples/$1/*.trimmed.unpaired.fq.gz
rm /mnt/gluster/bontrager2/popgen_samples/$1/*.fastq.gz
