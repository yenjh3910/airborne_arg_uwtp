#!/bin/bash

mkdir ~/contigs_kraken2

for i in $(<~/sample_list/clean_read_list)
do

	# kraken2
	kraken2 --db ~/db/kraken_db --threads 16 --report ~/contigs_kraken2/${i}_contigs_kraken2.report --output ~/contigs_kraken2/${i}_contigs_kraken2.output ~/megahit/megahit_individual/${i}.contigs.fa --use-names

done
