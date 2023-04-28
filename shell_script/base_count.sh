#!/bin/bash

for i in $(<~/sample_list/clean_read_list)
do
	echo "${i}_1"
	zcat ~/clean_read/${i}_1.fastq.gz | awk '{if(NR%4==2) print}' | tr -d "\n" | wc -c
	echo "${i}_2"
	zcat ~/clean_read/${i}_2.fastq.gz | awk '{if(NR%4==2) print}' | tr -d "\n" | wc -c
done