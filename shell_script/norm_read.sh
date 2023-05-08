#!/bin/bash

cd ~/clean_read

for i in $(<~/sample_list/clean_read_list)
do
	bbnorm.sh in1=${i}_1.fastq in2=${i}_2.fastq out1=n${i}_1.fastq out2=n${i}_2.fastq target=70
done