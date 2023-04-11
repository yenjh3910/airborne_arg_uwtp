#!/bin/bash

mkdir ~/contigs_diamond
mkdir ~/contigs_diamond/SARG
mkdir ~/contigs_diamond/MGE
mkdir ~/contigs_diamond/MRG

for i in $(<~/sample_list/clean_read_list)
do
	# Diamond
	# SARG blast
	diamond blastx -d ~/db/args_oap_db/SARG.dmnd -q ~/contigs_prodigal/${i}_contigs.nucl --id 70 -p 16 -e 1e-10 -k 1 --query-cover 70 -o  ~/contigs_diamond/SARG/${i}_contigs.SARG.dmnd

	# MGE blast
	diamond blastx -d  ~/db/MGE_db/MGE.dmnd -q ~/contigs_prodigal/${i}_contigs.nucl --id 70 -p 16 -e 1e-10 -k 1 --query-cover 70 -o ~/contigs_diamond/MGE/${i}_contigs.MGE.dmnd

	# MRG blast
	diamond blastx -d ~/db/BacMet_db/BacMet.dmnd -q ~/contigs_prodigal/${i}_contigs.nucl --id 70 -p 16 -e 1e-10 -k 1 --query-cover 70 -o ~/contigs_diamond/MRG/${i}_contigs.MRG.dmnd
done