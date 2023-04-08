#!/bin/bash

# Join the output files
mkdir ~/humann3/sample_join
humann_join_tables --input ~/humann3 --output ~/humann3/sample_join/genefamilies.tsv --file_name genefamilies
humann_join_tables --input ~/humann3 --output ~/humann3/sample_join/pathabundance.tsv --file_name pathabundance
humann_join_tables --input ~/humann3 --output ~/humann3/sample_join/pathcoverage.tsv --file_name pathcoverage

# gene family
## uniref90
mkdir ~/humann3/uniref90
humann_rename_table --input ~/humann/sample_join/genefamilies.tsv --names uniref90 --output ~/humann/uniref90/genefamilies_rename_uniref90.tsv
humann_renorm_table --input ~/humann/uniref90/genefamilies_rename_uniref90.tsv --output ~/humann/uniref90/genefamilies_rename_uniref90_cpm.tsv  --units cpm
humann_renorm_table --input ~/humann/uniref90/genefamilies_rename_uniref90.tsv --output ~/humann/uniref90/genefamilies_rename_uniref90_relab.tsv  --units relab

## Regroup to pfam
mkdir ~/humann3/pfam
humann_regroup_table -i ~/humann/sample_join/genefamilies.tsv -g uniref90_pfam -o ~/humann/pfam/genefamilies_regroup_pfam.tsv
humann_rename_table --input ~/humann/pfam/genefamilies_regroup_pfam.tsv  --names pfam --output ~/humann/pfam/genefamilies_rename_pfam.tsv
humann_renorm_table --input ~/humann/pfam/genefamilies_rename_pfam.tsv --output ~/humann/pfam/genefamilies_rename_pfam_cpm.tsv  --units cpm -s n

## Regroup to GO
mkdir ~/humann3/GO
humann_regroup_table -i ~/humann/sample_join/genefamilies.tsv -g uniref90_go -o ~/humann/GO/genefamilies_regroup_go.tsv
humann_rename_table --input ~/humann/GO/genefamilies_regroup_go.tsv  --names go --output ~/humann/GO/genefamilies_rename_go.tsv
humann_renorm_table --input ~/humann/GO/genefamilies_rename_go.tsv --output ~/humann/GO/genefamilies_rename_go_cpm.tsv  --units cpm -s n

## Regroup to EGGNOG
mkdir ~/humann3/EGGNOG
humann_regroup_table -i ~/humann/sample_join/genefamilies.tsv -g uniref90_eggnog -o ~/humann/EGGNOG/genefamilies_regroup_eggnog.tsv
humann_rename_table --input ~/humann/EGGNOG/genefamilies_regroup_eggnog.tsv  --names eggnog --output ~/humann/EGGNOG/genefamilies_rename_eggnog.tsv
humann_renorm_table --input ~/humann/EGGNOG/genefamilies_rename_eggnog.tsv --output ~/humann/EGGNOG/genefamilies_rename_eggnog_cpm.tsv  --units cpm -s n

## Regroup to KEGG
mkdir ~/humann3/KEGG
humann_regroup_table -i ~/humann/sample_join/genefamilies.tsv -g uniref90_ko -o ~/humann/KEGG/genefamilies_regroup_ko.tsv
humann_rename_table --input ~/humann/KEGG/genefamilies_regroup_ko.tsv  --names kegg-orthology --output ~/humann/KEGG/genefamilies_rename_kegg-orthology.tsv
humann_renorm_table --input ~/humann/KEGG/genefamilies_rename_kegg-orthology.tsv --output ~/humann/KEGG/genefamilies_rename_kegg-orthology_cpm.tsv  --units cpm -s n


# path abundance
mkdir ~/humann3/path_abundance
humann_renorm_table --input ~/humann/sample_join/pathabundance.tsv --output ~/humann/path_abundance/pathabundance_cpm.tsv  --units cpm
humann_renorm_table --input ~/humann/sample_join/pathabundance.tsv --output ~/humann/path_abundance/pathabundance_relab.tsv  --units relab