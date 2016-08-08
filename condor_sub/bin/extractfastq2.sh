#!/bin/bash
grep -v ^@ /mnt/gluster/bontrager2/bwa/MAE_meta.sam | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > /mnt/gluster/bontrager2/bwa/MAE_meta_2.fastq 

