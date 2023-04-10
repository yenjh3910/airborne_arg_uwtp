#!/bin/bash

mkdir ~/contigs_prodigal

for i in $(<~/sample_list/clean_read_list)
do

	# Prodigal
	prodigal -i ~/megahit/megahit_individual/${i}.contigs.fa -o ~/contigs_prodigal/${i}_contigs.gff -a ~/contigs_prodigal/${i}_contigs.prot -d ~/contigs_prodigal/${i}_contigs.nucl -c -p meta -f gff

done
