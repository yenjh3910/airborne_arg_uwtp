#!/bin/bash

for i in $(<~/sample_list/clean_read_list)
do
	megahit -t 16 -m 0.99 -1 ~/clean_read/${i}_1.fastq -2 ~/clean_read/${i}_2.fastq \
    --min-contig-len 1000 -o ~/megahit/megahit_individual/${i}_contigs --presets meta-large
    mv  ~/megahit/megahit_individual/${i}_contigs/final.contigs.fa ~/megahit/megahit_individual/${i}.contigs.fa
done