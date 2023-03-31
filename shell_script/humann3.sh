#!/bin/bash

for i in $(<~/sample_list/clean_read_list)
do

cat ~/clean_read/${i}_1.fastq.gz  ~/clean_read/${i}_2.fastq.gz > ~/clean_read/${i}.fastq.gz

# Run humann3
humann -i ~/clean_read/${i}.fastq.gz -o ~/humann3 --input-format fastq.gz --threads 16 --metaphlan-options "--bowtie2db db/humann3_db/metaphlan_database_humann_compatible --index mpa_vJan21_CHOCOPhlAnSGB_202103"

rm ~/clean_read/${i}.fastq.gz

done