#!/bin/bash
bwa/bwa mem -M -t 20  /mnt/gluster/bontrager2/gen_consensus_test/green /mnt/gluster/bontrager2/CBE_R1.fastq.gz /mnt/gluster/bontrager2/CBE_R2.fastq.gz > /mnt/gluster/bontrager2/bwa/BWAgreen.test.sam

