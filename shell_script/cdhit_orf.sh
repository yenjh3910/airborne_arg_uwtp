#!/bin/bash

mkdir ~/contigs_cdhit

for i in $(<~/sample_list/clean_read_list)
do

	# CD-HIT
	cd-hit -i ~/contigs_prodigal/${i}_contigs.nucl -o ~/contigs_cdhit/${i}.nucl.uniq -M 0 -T 0 -l 250 -s 0.9 -c 0.9
	cd-hit -i ~/contigs_prodigal/${i}_contigs.prot -o ~/contigs_cdhit/${i}.prot.uniq -M 0 -T 0 -l 80 -c 0.9

done
