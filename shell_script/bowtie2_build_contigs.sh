#!/bin/bash

mkdir ~/contigs_bowtie2
mkdir ~/contigs_bowtie2/SARG
mkdir ~/contigs_bowtie2/MGE
mkdir ~/contigs_bowtie2/VF
mkdir ~/contigs_bowtie2/SARG/index
mkdir ~/contigs_bowtie2/MGE/index
mkdir ~/contigs_bowtie2/VF/index

for i in $(<~/sample_list/clean_read_list)
do
	#Build index
	bowtie2-build ~/extract_contigs/SARG/${i}_SARG.fa ~/contigs_bowtie2/SARG/index/${i}_SARG.index
	bowtie2-build ~/extract_contigs/MGE/${i}_MGE.fa ~/contigs_bowtie2/MGE/index/${i}_MGE.index
	bowtie2-build ~/extract_contigs/VF/${i}_VF.fa ~/contigs_bowtie2/VF/index/${i}_VF.index
done
