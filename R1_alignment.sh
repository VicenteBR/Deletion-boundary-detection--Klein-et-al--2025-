#!/bin/bash
for i in *_1.fq.gz;
do

read1=$(echo -e $i | sed 's/_1.fq.gz//g')
outfile=$(echo -e $read1 | sed 's/_.*//g')

#Alignment to E.coli BL21AI genome - R1 with splice
        echo -e "Processing ${Outfile} - R1 with splice. \n"
	date +"%T"
                        hisat2 -x /mnt/e/nk_deletion/R1_R2_splice_single/index/ecomplete \
                        -U $i \
                        -S /mnt/e/nk_deletion/R1_R2_splice_single/R1/"${read1}_R1.sam" \
                        --no-softclip \
                        -p 60 \


done
