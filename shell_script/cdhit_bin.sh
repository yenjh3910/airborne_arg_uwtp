#!/bin/bash

mkdir ~/bins_cdhit

for i in {1..112}

do
	# CD-HIT
	cd-hit -i ~/bins_prodigal/bins.${i}.nucl -o ~/bins_cdhit/bin.${i}.nucl.uniq -M 0 -T 0 -l 250 -s 0.9 -c 0.9
	cd-hit -i ~/bins_prodigal/bins.${i}.prot -o ~/bins_cdhit/bin.${i}.prot.prot -M 0 -T 0 -l 80 -c 0.9
done
