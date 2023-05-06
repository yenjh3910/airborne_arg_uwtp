#!/bin/bash

mkdir ~/contigs_diamond
mkdir ~/contigs_diamond/SARG
mkdir ~/contigs_diamond/MGE
mkdir ~/contigs_diamond/MRG
mkdir ~/contigs_diamond/VF

for i in $(<~/sample_list/clean_read_list)
do
	# Diamond
	# SARG blast
	diamond blastx -d ~/db/args_oap_db/SARG.dmnd -q ~/contigs_cdhit/${i}_contigs.nucl.uniq --id 50 -p 16 -e 1e-10 -k 1 --query-cover 50 -o  ~/contigs_diamond/SARG/${i}_contigs.SARG.dmnd

	# MGE blast
	diamond blastx -d  ~/db/MGE_db/MGE.dmnd -q ~/contigs_cdhit/${i}_contigs.nucl.uniq --id 50 -p 16 -e 1e-10 -k 1 --query-cover 50 -o ~/contigs_diamond/MGE/${i}_contigs.MGE.dmnd

	# MRG blast
	diamond blastx -d ~/db/BacMet_db/BacMet.dmnd -q ~/contigs_cdhit/${i}_contigs.nucl.uniq --id 50 -p 16 -e 1e-10 -k 1 --query-cover 50 -o ~/contigs_diamond/MRG/${i}_contigs.MRG.dmnd

	# VFDB blast
	diamond blastx -d ~/db/vfdb/VFDB.dmnd -q ~/contigs_cdhit/${i}_contigs.nucl.uniq --id 50 -p 16 -e 1e-10 -k 1 --query-cover 50 -o ~/contigs_diamond/VF/${i}_contigs.VF.dmnd
done
