#!/bin/bash
mkdir -p bam_outs
mkdir -p bam_outs/standard_bam

for i in *.sam;
do
        outfile=$(echo $i | sed 's/_.*//g')
        echo -e "Processing ${outfile}. \n"
        echo -e "Converting and sorting to standard bam file. \n"

                samtools view -Sbq 2 -@ 24 $i -o t1.bam
                samtools sort -@ 24 t1.bam -o bam_outs/standard_bam/$outfile"_sorted.bam"
                samtools index -@ 24 bam_outs/standard_bam/$outfile"_sorted.bam"

done
