#!/bin/bash

mkdir ~/assembly_individual_coverage
mkdir ~/assembly_individual_coverage/index
mkdir ~/assembly_individual_coverage/sam
mkdir ~/assembly_individual_coverage/map

for i in $(<~/sample_list/clean_read_list)
do
	# Create bowtie2 index database
	bowtie2-build ~/megahit/megahit_individual/${i}.contigs.fa \
  	~/assembly_individual_coverage/index/${i}.contigs.index
	# Map reads
	bowtie2 -x ~/assembly_individual_coverage/index/${i}.contigs.index \
	-1 ~/clean_read/${i}_1.fastq.gz -2 ~/clean_read/${i}_2.fastq.gz \
	-S ~/assembly_individual_coverage/sam/${i}.contigs.sam  -p 16
	# Print mapping read number and generate mapping file
	pileup.sh in=~/assembly_individual_coverage/sam//${i}.contigs.sam \
	out=~/assembly_individual_coverage/map/${i}.contigs.sam.map.txt

	rm ~/assembly_individual_coverage/sam/${i}.contigs.sam
done
