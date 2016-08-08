#!/bin/bash
SPAdes-3.5.0-Linux/bin/spades.py --dataset /mnt/gluster/bontrager2/condor_output/2015-09-08_MME/corrected/corrected.yaml \
   --only-assembler -k 33,55,77 -t 20 -m 150 -o /mnt/gluster/bontrager2/condor_output/2015-09-08_MME/assembly
