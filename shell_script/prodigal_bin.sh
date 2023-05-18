#!/bin/bash

mkdir ~/bins_prodigal

for i in {1..112}

do
	# Prodigal
	prodigal -i ~/metawrap_run/bin_refinement/metawrap_50_10_bins/bin.${i}.fa -o ~/bins_prodigal/bins.${i}.prodigalout -a ~/bins_prodigal/bins.${i}.prot -d ~/bins_prodigal/bins.${i}.nucl -c -p meta
done
