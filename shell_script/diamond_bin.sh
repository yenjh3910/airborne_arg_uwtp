#!/bin/bash

mkdir ~/bins_diamond
mkdir ~/bins_diamond/SARG
mkdir ~/bins_diamond/MGE
mkdir ~/bins_diamond/MRG
mkdir ~/bins_diamond/VF


for i in {1..112}
do
	# Diamond
	# SARG blast
	diamond blastx -d ~/db/args_oap_db/SARG.dmnd -q ~/bins_cdhit/bin.${i}.nucl.uniq --id 50 -p 16 -e 1e-10 -k 1 --query-cover 50 -o  ~/bins_diamond/SARG/bin.${i}.SARG.dmnd

	# MGE blast
	diamond blastx -d  ~/db/MGE_db/MGE.dmnd -q ~/bins_cdhit/bin.${i}.nucl.uniq --id 50 -p 16 -e 1e-10 -k 1 --query-cover 50 -o ~/bins_diamond/MGE/bin.${i}.MGE.dmnd

	# MRG blast
	diamond blastx -d ~/db/BacMet_db/BacMet.dmnd -q ~/bins_cdhit/bin.${i}.nucl.uniq --id 50 -p 16 -e 1e-10 -k 1 --query-cover 50 -o ~/bins_diamond/MRG/bin.${i}.MRG.dmnd

	# VFDB blast
	diamond blastx -d ~/db/vfdb/VFDB.dmnd -q ~/bins_cdhit/bin.${i}.nucl.uniq --id 50 -p 16 -e 1e-10 -k 1 --query-cover 50 -o ~/bins_diamond/VF/bin.${i}.VF.dmnd
done
