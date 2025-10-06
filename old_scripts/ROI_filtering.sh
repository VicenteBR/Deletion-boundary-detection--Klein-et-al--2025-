#!/bin/bash

for i in *.bam;
do
    #Select reads that are found inside the region of interest (300000 to 360000) in the E. coli BL21Ai genome
    samtools view -h $i | awk 'BEGIN {OFS="\t"} /^@/ || ($3 == "CP047231.1" && $4 >= 300000 && $4 + length($10) - 1 <= 360000)' > t1.sam
    samtools view -Sb t1.sam > down.bam
    samtools sort down.bam > filtered_bams/"filtered_"$i
    samtools index filtered_bams/"filtered_"$i
    rm t1.bam down.bam t1.sam

done
