#!/bin/bash
SPAdes-3.5.0-Linux/bin/spades.py --pe1-1 /mnt/gluster/bontrager2/bwa/MME_meta_1.fastq \
--pe1-2 /mnt/gluster/bontrager2/bwa/MME_meta_2.fastq -o 2015-09-08_MME -t 25 \
  -m 200 --only-error-correction
cp -R 2015-09-08_MME /mnt/gluster/bontrager2/condor_output/
