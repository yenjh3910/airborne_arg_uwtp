#!/bin/bash

mkdir ~/prokka_annotation

for i in {1..112}
do
	# SARG blast
	prokka ~/metawrap_run/bin_refinement/metawrap_50_10_bins/bin.${i}.fa --outdir ~/prokka_annotation/bin.${i} -prefix bin.${i}

done
