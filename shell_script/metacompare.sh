#!/bin/bash

cd ~/MetaCompare

for i in $(<~/sample_list/clean_read_list)
do
    # move contgs file to ~/metacompare
    mv ~/megahit/megahit_individual/${i}.contigs.fa ~/MetaCompare 
    # move prodigal nucl to ~/metacompare
    mv ~/contigs_prodigal/${i}_contigs.nucl ~/MetaCompare
    
    # run metacompare
    ~/MetaCompare/metacmp.py -c ${i}.contigs.fa -g ${i}_contigs.nucl -t 16
    
    #move file back
    mv ~/MetaCompare/${i}.contigs.fa ~/megahit/megahit_individual
    mv ~/MetaCompare/${i}_contigs.nucl  ~/contigs_prodigal
done