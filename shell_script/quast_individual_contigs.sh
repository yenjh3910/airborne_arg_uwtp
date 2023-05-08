#!/bin/bash

for i in $(<~/sample_list/clean_read_list)
do
    mkdir ~/megahit/megahit_individual/${i}_contigs_quast
    ~/quast/quast.py ~/megahit/megahit_individual/${i}.contigs.fa -o  ~/megahit/megahit_individual/${i}_contigs_quast
done