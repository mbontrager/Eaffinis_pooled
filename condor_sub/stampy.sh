#!/bin/bash
# Add the /mnt/gluster software location to the job's PATH
export PATH=/mnt/gluster/bontrager2/bin/Python2.7/bin:$PATH

# Run stampy
python2.7 /mnt/gluster/bontrager2/bin/stampy-1.0.28/stampy.py \
-g /mnt/gluster/bontrager2/eaff_masked/eaff_masked \
-h /mnt/gluster/bontrager2/eaff_masked/eaff_masked \
-t 32 \
--substitutionrate=0.07 \
--insertsize=450 \
--bamkeepgoodreads \
-M /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.bam \
-o /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.stampy.sam

#Process
#rm /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.bam
samtools/samtools view -bS /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.stampy.sam \
-o /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.stampy.bam
rm /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.stampy.sam
samtools/samtools sort /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.stampy.bam \
/mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.stampy.sorted
rm /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.stampy.bam
samtools/samtools index /mnt/gluster/bontrager2/popgen_samples/$1/$1.bwa.stampy.sorted.bam
