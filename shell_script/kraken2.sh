#!/bin/bash

for i in $(<~/sample_list/clean_read_list)
do

	# kraken2
	kraken2 --db ~/db/kraken_db --threads 16 --report ~/kraken2/${i}_kraken2.report --output ~/kraken2/${i}_kraken2.output --paired ~/clean_read/${i}_1.fastq.gz ~/clean_read/${i}_2.fastq.gz
	
	# bracken
	bracken -d ~/db/kraken_db -i ~/kraken2/${i}_kraken2.report -o ~/kraken2/${i}.D.bracken -w ~/kraken2/${i}.D.bracken.report -r 150 -l D
	bracken -d ~/db/kraken_db -i ~/kraken2/${i}_kraken2.report -o ~/kraken2/${i}.P.bracken -w ~/kraken2/${i}.P.bracken.report -r 150 -l P
	bracken -d ~/db/kraken_db -i ~/kraken2/${i}_kraken2.report -o ~/kraken2/${i}.C.bracken -w ~/kraken2/${i}.C.bracken.report -r 150 -l C
	bracken -d ~/db/kraken_db -i ~/kraken2/${i}_kraken2.report -o ~/kraken2/${i}.O.bracken -w ~/kraken2/${i}.O.bracken.report -r 150 -l O
	bracken -d ~/db/kraken_db -i ~/kraken2/${i}_kraken2.report -o ~/kraken2/${i}.F.bracken -w ~/kraken2/${i}.F.bracken.report -r 150 -l F
	bracken -d ~/db/kraken_db -i ~/kraken2/${i}_kraken2.report -o ~/kraken2/${i}.G.bracken -w ~/kraken2/${i}.G.bracken.report -r 150 -l G
	bracken -d ~/db/kraken_db -i ~/kraken2/${i}_kraken2.report -o ~/kraken2/${i}.S.bracken -w ~/kraken2/${i}.S.bracken.report -r 150 -l S

done
