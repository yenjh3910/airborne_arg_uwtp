#!/bin/bash

mkdir ~/contigs_bowtie2/SARG/sam
mkdir ~/contigs_bowtie2/MGE/sam
mkdir ~/contigs_bowtie2/VF/sam
mkdir ~/contigs_bowtie2/SARG/coverage
mkdir ~/contigs_bowtie2/MGE/coverage
mkdir ~/contigs_bowtie2/VF/coverage

for i in $(<~/sample_list/clean_read_list)
do
	# SARG
	## Bowtie2 mapping
	bowtie2 -x ~/contigs_bowtie2/SARG/index/${i}_SARG.index -1 ~/clean_read/${i}_1.fastq.gz -2 ~/clean_read/${i}_2.fastq.gz -S ~/contigs_bowtie2/SARG/sam/${i}_SARG.sam -p 16
	## Coverage calculation
	pileup.sh in=~/contigs_bowtie2/SARG/sam/${i}_SARG.sam out=~/contigs_bowtie2/SARG/coverage/${i}_SARG.sam.map.txt
	rm ~/contigs_bowtie2/SARG/sam/${i}_SARG.sam

	# MGE
	## Bowtie2 mapping
	bowtie2 -x ~/contigs_bowtie2/MGE/index/${i}_MGE.index -1 ~/clean_read/${i}_1.fastq.gz -2 ~/clean_read/${i}_2.fastq.gz -S ~/contigs_bowtie2/MGE/sam/${i}_MGE.sam -p 16
	## Coverage calculation
	pileup.sh in=~/contigs_bowtie2/MGE/sam/${i}_MGE.sam out=~/contigs_bowtie2/MGE/coverage/${i}_MGE.sam.map.txt
	rm ~/contigs_bowtie2/MGE/sam/${i}_MGE.sam

	# VF
	## Bowtie2 mapping
	bowtie2 -x ~/contigs_bowtie2/VF/index/${i}_VF.index -1 ~/clean_read/${i}_1.fastq.gz -2 ~/clean_read/${i}_2.fastq.gz -S ~/contigs_bowtie2/VF/sam/${i}_VF.sam -p 16
	## Coverage calculation
	pileup.sh in=~/contigs_bowtie2/VF/sam/${i}_VF.sam out=~/contigs_bowtie2/VF/coverage/${i}_VF.sam.map.txt
	rm ~/contigs_bowtie2/VF/sam/${i}_VF.sam
done