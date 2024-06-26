# Airborne ARGs in UWTP
A comprehensive metagenomic pipeline for airborne ARGs (Antibiotic Resistance Genes) in UWTP (Urban wastewater Treatment Plant)

# Table of Contents
- [Environment Setup](https://github.com/yenjh3910/airborne_arg_uwtp#environment-setup)
    - [Anaconda Installation](https://github.com/yenjh3910/airborne_arg_uwtp#anaconda-installation)
- [Reference-based Analysis](https://github.com/yenjh3910/airborne_arg_uwtp#reference-based-analysis)
    - [Taxanomic Profile](https://github.com/yenjh3910/airborne_arg_uwtp#taxanomic-profile)
    - [ARGs Profile](https://github.com/yenjh3910/airborne_arg_uwtp#args-profile)
    - [Functional Profile](https://github.com/yenjh3910/airborne_arg_uwtp#functional-profile)
- [Assembly-based Analysis](https://github.com/yenjh3910/airborne_arg_uwtp#assembly-based-analysis)
    - [Contigs Assembly](https://github.com/yenjh3910/airborne_arg_uwtp#contigs-assembly)
    - [Taxanomic Assignment of Assembly Contigs](https://github.com/yenjh3910/airborne_arg_uwtp#taxanomic-assignment-of-assembly-contigs)
    - [Gene Alignment to Assembly Contigs](https://github.com/yenjh3910/airborne_arg_uwtp#gene-alignment-to-assembly-contigs)
    - [Calculate coverage of aligning contigs](https://github.com/yenjh3910/airborne_arg_uwtp#calculate-coverage-of-aligning-contigs)
    - [Prediction of Plasmid Sequences](https://github.com/yenjh3910/airborne_arg_uwtp#prediction-of-plasmid-sequences)
    - [Risk Assessment](https://github.com/yenjh3910/airborne_arg_uwtp#risk-assessment)
- [Binning-based Analysis](https://github.com/yenjh3910/airborne_arg_uwtp#binning-based-analysis)
    - [Contigs Co-assembly](https://github.com/yenjh3910/airborne_arg_uwtp/tree/master#contigs-co-assembly)
    - [Binning](https://github.com/yenjh3910/airborne_arg_uwtp#binning)
    - [Taxonomy Classification of Bin](https://github.com/yenjh3910/airborne_arg_uwtp#taxonomy-classification-of-bin)
    - [Gene Alignment to Bin](https://github.com/yenjh3910/airborne_arg_uwtp#gene-alignment-to-bin)
# Environment Setup
## [Anaconda Installation](https://www.anaconda.com/products/distribution)
```
$ sudo apt update && upgrade
$ sudo apt install curl
$ mkdir tmp
$ cd tmp

# Install
$ curl -O https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh
$ bash Anaconda3-2022.10-Linux-x86_64.sh
$ vi ~/.bash_profile     # Add 'export PATH=$PATH:/home/yen/anaconda3/bin' in last row
$ source ~/.bashrc
$ conda info

$ cd
$ rm -r tmp
```
### [Bioconda Installation](https://bioconda.github.io/)
```
# Bioconda set up channels
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
$ conda config --set channel_priority strict
# View channel
$ conda config --show channels
```
# Reference-based Analysis
## Taxanomic Profile
### [Kraken2](https://github.com/DerrickWood/kraken2/wiki/Manual)
#### Create environment and install
```
$ conda create --name kraken2
$ conda activate kraken2

# Install
$ conda install kraken2

# Version
$ kraken2 --version
$ conda update kraken2 -c bioconda
```
#### Download NCBI taxonomic information
```
$ cd ~/db
$ kraken2-build --download-taxonomy --threads 16 --db kraken_db
```
#### Debug  
1. Refer to https://qiita-com.translate.goog/kohei-108/items/ce5fdf10c11d1e7ca15b?_x_tr_sl=ja&_x_tr_tl=zh-TW&_x_tr_hl=zh-TW&_x_tr_pto=sc :  
Change rsync_from_ncbi.pl from https://translate.google.com/website?sl=ja&tl=zh-TW&hl=zh-TW&prev=search&u=https://github.com/DerrickWood/kraken2/issues/465%23issuecomment-870726096  

2. https://github.com/DerrickWood/kraken2/issues/518  
rsync_from_ncbi.pl  Line 46:  
```
if (! ($full_path =~ s#^ftp://${qm_server}${qm_server_path}/##))  to 
if (! ($full_path =~ s#^https://${qm_server}${qm_server_path}/##))
```
#### Download partial library
```
$ kraken2-build --download-library bacteria --threads 16 --db kraken_db
```
#### Create DB from library
```
$ cd ~/db/kraken_db
$ kraken2-build --build --threads 16 --db ./ --max-db-size 52000000000
```
### [Bracken](https://github.com/jenniferlu717/Bracken)
#### Installation
```
# Install prerequisite
$ sudo apt-get install build-essential

# Install Bracken
$ git clone https://github.com/jenniferlu717/Bracken.git
$ cd Bracken
$ bash install_bracken.sh
$ cd src/ && make

# Add dictionary to path
$ nano ~/.bashrc
# Add to last row
$ export PATH="~/Bracken:$PATH"
$ export PATH="~/Bracken/src:$PATH"
# Load the new $PATH
$ source ~/.bashrc

# Test 
$ bracken -h
```
#### Build bracken library
```
$ bracken-build -d ~/db/kraken_db -t 16 -k 35 -l 150
```
#### Alternative (download pre-built database)
https://benlangmead.github.io/aws-indexes/k2

#### Usage 
 ```
$ ~/shell_script/kraken2.sh
 ```
### R script after Taxanomic Profile
```
$ cd ./taxonomy
$ Rscript phylum_barplot.R
$ Rscript clss_barplot.R
$ Rscript order_barplot.R
$ Rscript family_barplot.R
$ Rscript genus_barplot.R
$ Rscript species_barplot.R
$ Rscript taxa_alpha_diversity.R
```
## ARGs Profile
### [ARGs-OAP](https://github.com/xinehc/args_oap)
#### Installation
```
$ conda create -n args_oap
$ conda activate args_oap
$ conda config --set channel_priority flexible
$ conda install -c bioconda -c conda-forge args_oap
```
#### ARGs analysis
```
# Stage one
$ args_oap stage_one -i ~/clean_read -o ~/args_oap/ARG/stage_one_output -f fastq -t 16

# Stage two (e_value: 1e-7, identity: 80, aa_length,25)
$ args_oap stage_two -i ~/args_oap/ARG/stage_one_output -o ~/args_oap/ARG/stage_two_output -t 16
```
#### MGEs (Mobile Genetic Elements) analysis
Database: https://github.com/KatariinaParnanen/MobileGeneticElementDatabase
1. Covert nucleotide acid to amino acid under fasta format
2. Create MGE_structure.txt manually or use curated structure already made
Curated MGE database can be found [here](https://github.com/yenjh3910/airborne_arg_uwtp/blob/master/MGE/MGE_structure/MGE_curated_structure.txt)
```
$ mkdir ~/args_oap/MGE
$ cd ~/args_oap/MGE

### Amino acid
# (Optional)
$ echo '>level1' | cat - MGEs_database_aa.fasta | grep '^>' | cut -d ' ' -f 1 | cut -c2- > MGE_AA_structure.txt

# The database should be indexed manually (protein or nucleotide, in fasta)
$ args_oap make_db -i MGEs_database_aa.fasta

# Stage one
$ args_oap stage_one -i ~/clean_read -o ~/args_oap/MGE/AA_stage_one_output -f fastq -t 16 --database ~/args_oap/MGE/MGEs_database_aa.fasta

# Stage two
$ args_oap stage_two -i ~/args_oap/MGE/AA_stage_one_output -o ~/args_oap/MGE/AA_stage_two_output -t 16 --database ~/args_oap/MGE/MGEs_database_aa.fasta --structure1 ~/args_oap/MGE/MGE_curated_structure.txt

### DNA blast
# The database should be indexed manually (protein or nucleotide, in fasta)
$ args_oap make_db -i MGEs_database_dna.fasta

# Stage one
$ args_oap stage_one -i ~/clean_read -o ~/args_oap/MGE/DNA_stage_one_output -f fastq -t 16 --database ~/args_oap/MGE/MGEs_database_dna.fasta

# Stage two
$ args_oap stage_two -i ~/args_oap/MGE/DNA_stage_one_output -o ~/args_oap/MGE/DNA_stage_two_output -t 16 --database ~/args_oap/MGE/MGEs_database_dna.fasta --structure1 ~/args_oap/MGE/MGE_curated_structure.txt
```
### MRGs (Metal Resistance Genes) analysis
Database: http://bacmet.biomedicine.gu.se/  
Curated metal structure file can be found [here](https://github.com/yenjh3910/airborne_arg_uwtp/blob/master/BacMet/BacMet_structure/metal_only_structure.txt)
```
# The database should be indexed manually (protein or nucleotide, in fasta)
$ args_oap make_db -i BacMet_exp_metal.fasta

# Stage one
$ args_oap stage_one -i ~/clean_read -o ~/args_oap/BacMet/stage_one_output -f fastq -t 16 --database ~/args_oap/BacMet/BacMet_exp_metal.fasta

# Stage two
$ args_oap stage_two -i ~/args_oap/BacMet/stage_one_output -o ~/args_oap/BacMet/stage_two_output -t 16 --database ~/args_oap/BacMet/BacMet_exp_metal.fasta --structure1 ~/args_oap/BacMet/metal_only_structure.txt

# Since default parameter in stage two is too strict for MGE, following parameters (--e 1e-5 --id 70) were used: 
$ args_oap stage_two -i ~/args_oap/BacMet/stage_one_output -o ~/args_oap/BacMet/stage_two_output_evalue-5_id70 --e 1e-5 --id 70 -t 16 --database ~/args_oap/BacMet/BacMet_exp_metal.fasta --structure1 ~/args_oap/BacMet/metal_only_structure.txt
```
### R script after ARGs Profile
```
$ cd ./ARG
$ Rscript ARG_type_cell_abundance.R
$ Rscript ARG_type_16S_abundance.R
$ Rscript ARG_type_ppm_abundance.R
$ Rscript ARG_type_rpkm_abundance.R
$ Rscript ARG_type_tpm_abundance.R
$ Rscript ARG_subtype_cell_abundance.R
$ Rscript ARG_type_circos.R
$ Rscript ARG_subtype_venn_diagram.R
$ Rscript ARG_UpSet.R
$ Rscript ARG_mechanism.R
$ Rscript ARG_subtype_cell_abundance.R
cd ../MGE
$ Rscript MGE_type_cell_abundance.R
cd ../BacMet
$ Rscript MRG_type_cell_abundance.R
cd ../taxonomy
$ Rscript ARG_species_procrustes.R
cd ../reference_based_correlation
$ Rscript Correlation_ARG_MGE.R
$ Rscript Correlation_MRG_MGE.R
cd ../LEFSe
$ Rscript arg_lefse.R
```
## Functional Profile
### [HUMAnN 3.0](https://github.com/biobakery/humann/tree/8d69f3c84ca7bfd7519ced7fcf94b8356c915090)
#### Installation
```
# Create environment
$ conda create --name humann3
$ conda activate humann3

# Add channel
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
$ conda config --add channels biobakery

# Install
$ conda install -c conda-forge -c bioconda -c biobakery metaphlan=4.0.3
$ conda install -c biobakery humann

# Test
$ cd ~/clean_read
$ metaphlan <sample_1.fastq.gz,sample_2.fastq.gz> --bowtie2out sample_metaphlan.bowtie2.bz2 --input_type fastq --nproc 16 > sample_metaphlan.txt
$ humann_test
```
#### Download database
```
$ humann_databases --download chocophlan full ~/db/humann3_db --update-config yes # ChocoPhlAn database
$ humann_databases --download uniref uniref90_diamond ~/db/humann3_db --update-config yes # Translated search databases
$ humann_databases --download utility_mapping full  ~/db/humann3_db --update-config yes # HUMAnN 3.0 utility mapping files
```
#### Update path
```
$ humann_config --update database_folders nucleotide ~/db/humann3_db/chocophlan
$ humann_config --update database_folders protein ~/db/humann3_db/uniref
$ humann_config --update database_folders utility_mapping ~/db/humann3_db/utility_mapping
```
#### [Bug when running humann3](https://forum.biobakery.org/t/metaphlan-v4-0-2-and-huma-3-6-metaphlan-taxonomic-profile-provided-was-not-generated-with-the-expected-database/4296/22): 
```
config.metaphlan_v3_db_version+" or "+metaphlan_v4_db_version+" . Please update your version of MetaPhlAn to at least v3.0."
NameError: name 'metaphlan_v4_db_version' is not defined

# Solve:
## Install the database in a folder outside the Conda environment with a specific version.
$ metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db db/humann3_db/metaphlan_database_humann_compatible
```
#### Usage
```
# Run
$ ~/shell_script/humann3.sh

# Post processing
$ ~/shell_script/humann3_post_processing.sh
```
### R script after Functional Profile
```
$ Rscript humann_aerOnly_heatmap.R
$ Rscript humann3_heatmap_without_stress.R
```
# Assembly-based Analysis
## Contigs Assembly
### [MEGAHIT](https://github.com/voutcn/megahit) & [Quast](https://github.com/ablab/quast)
```
# Create environment
$ conda create -n megahit
$ conda activate megahit

# Installation
## Megahit
$ mamba install -c bioconda megahit
## Quast
$ git clone https://github.com/ablab/quast.git
$ cd ~/quast
$ ./setup.py install_full
$ sudo apt-get update && sudo apt-get install -y pkg-config libfreetype6-dev libpng-dev python3-matplotlib

# Run
## Contigs assembly
$ ~/shell_script/assmebly_individual.sh
## Quality check
$ ~/shell_script/quast_individual_contigs.sh

# View mapping read to contigs
$ conda activate diamond
$ ~/shell_script/assembly_individual_coverage.sh
```
## Taxanomic Assignment of Assembly Contigs
Use kraken2 to assign the toxanomy to individual assmebly contigs  
#### Usage
```
$ conda activate kraken2
$ ~/shell_script/kraken2_contigs.sh
```
## Gene Alignment to Assembly Contigs
### [Prodigal](https://github.com/hyattpd/Prodigal) & [CD-HIT](https://github.com/weizhongli/cdhit)
Gene prediction & reductant sequence remove
```
# Create environment
$ conda create -n prodigal
$ conda activate prodigal

# Installation
$ conda install -c bioconda prodigal
$ conda install -c bioconda cd-hit

# ORF perdiction
$ ~/shell_script/prodigal_contigs.sh

# Remove reductant sequence
$ ~/shell_script/cdhit_orf.sh
```
### [Diamond](https://github.com/bbuchfink/diamond)
BLAST sequence
```
# Create environment
$ conda create -n diamond python=3.10
$ conda activate diamond

# Installation
$ conda install -c bioconda diamond

# Make database
$ diamond makedb --in ~/db/args_oap_db/sarg.fasta --db ~/db/args_oap_db/SARG.dmnd
$ diamond makedb --in ~/db/MGE_db/MGEs_database_aa.fasta --db ~/db/MGE_db/MGE.dmnd
$ diamond makedb --in ~/db/BacMet_db/BacMet_exp_metal.fasta --db ~/db/BacMet_db/BacMet.dmnd
$ diamond makedb --in ~/db/vfdb/VFDB_setB_pro.fa --db ~/db/vfdb/VFDB.dmnd

# Run
$ ~/shell_script/diamond_contigs.sh
```
### R script after Gene Alignment to Assembly Contigs
```
$ cd ./contigs

# Merge blast contigs with kraken2 taxonomy
$ Rscript contigs_kraken2.R

# QC & merge blast contigs with database
$ Rscript contigs_ARG_diamond.R
$ Rscript contigs_MGE_diamond.R
$ Rscript contigs_VF_diamond.R

# Extract contigs as mapping reference for coverage calculation
$ Rscript extract_SARG_contigs.R
$ Rscript extract_MGE_contigs.R
$ Rscript extract_VF_contigs.R
```
## Calculate coverage of aligning contigs
### [Bowtie2](https://github.com/BenLangmead/bowtie2) & [BBMap](https://github.com/BioInfoTools/BBMap)
```
# Enter environment
$ conda activate diamond

# Installation
$ conda install -c bioconda bowtie2
$ conda install -c bioconda bbmap

# Build bowtie2 index
$ ~/shell_script/bowtie2_build_contigs.sh

# Bowtie2 mapping & coverage calculation
$  ~/shell_script/coverage_contigs.sh

# Post processing (Coverage normalization & ORF annotation)
$ ~/shell_script/base_count.sh
Then create /airborne_arg_uwtp_result/read_base_count.xlsx manually
```
### R script after Calculate coverage of aligning contigs
```
$ cd ./contigs
$ Rscript ARG_coverage.R    # Plot ARG coverage in each ARB taxonomic level
$ Rscript MGE_coverage.R
$ Rscript VF_coverage.R
$ Rscript extract_ORF_position.R
$ Rscript bind_ARG_MGE_VF_coverage.R    # Plot coverage correlation of ARG-MGE & co-occurence of ARG-MGE in gggenes
$ Rscript ARG_coverage_sankey.R
```
## Prediction of Plasmid Sequences
### [PlasFlow](https://github.com/smaegol/PlasFlow)
```
```
## Risk Assessment
### [MetaCompare](https://github.com/minoh0201/MetaCompare)
```
# Installation
$ sudo apt-get update
$ sudo apt-get install python3-biopython
$ sudo apt-get install python3-pandas
$ git clone https://github.com/minoh0201/MetaCompare
$ cd MetaCompare
$ mkdir BlastDB
$ cd BlastDB
$ wget http://bench.cs.vt.edu/ftp/data/metacomp/BlastDB.tar.gz
$ tar -zxvf BlastDB.tar.gz
$ cd ..
$ ./metacmp.py

# Run
$ ~/shell_script/metacompare.sh
```
### R script after metacompare
```
$ cd ./metacompare
$ Rscript metacompare_plot.R
```
# Binning-based Analysis
## Contigs Co-assembly
### Reads normalization
```
# Enter environment
$ conda activate diamond

# Unzip fastq files
$ gunzip ~/clean_read/*_1.fastq.gz ~/clean_read/*_2.fastq.gz

# Run
$ ~/shell_script/norm_read.sh
```
### Run co-assembly
```
# Enter environment
$ conda activate megahit

# Concatenate reads
$ cat ~/clean_read/nA*_1.fastq > ~/clean_read/naer_reads_1.fastq
$ cat ~/clean_read/nA*_2.fastq > ~/clean_read/naer_reads_2.fastq

# Run co-assembly
$ megahit -t 16 -m 0.99 -1 ~/clean_read/naer_reads_1.fastq -2 ~/clean_read/naer_reads_2.fastq \
  --min-contig-len 1000 -o ~/megahit/megahit_coassembly/aeration_AT_ARP --presets meta-large

# Quality check
$ ~/quast/quast.py ~/megahit/megahit_coassembly/aeration_AT_ARP/final.contigs.fa \
  -o ~/megahit/megahit_coassembly/aeration_AT_ARP/coassembly_aeration_quast
```
### Get read and base number mapped to co-assembly contigs
```
# Enter environment
$ conda activate diamond

# Create bowtie2 index database
$ mkdir ~/megahit/megahit_coassembly/aeration_AT_ARP/breadth_coverage
$ bowtie2-build ~/megahit/megahit_coassembly/aeration_AT_ARP/final.contigs.fa \
  ~/megahit/megahit_coassembly/aeration_AT_ARP/breadth_coverage/coassembly_aer_contig.index

# AT
# Map reads
bowtie2 -x ~/megahit/megahit_coassembly/aeration_AT_ARP/breadth_coverage/coassembly_aer_contig.index \
  -1 ~/clean_read/AT*_1.fastq.gz -2 ~/clean_read/AT*_2.fastq.gz \
  -S ~/megahit/megahit_coassembly/aeration_AT_ARP/breadth_coverage/coassembly_AT_to_aer.sam.map  -p 16

# Print mapping read number and generate mapping file
$ pileup.sh in=~/megahit/megahit_coassembly/aeration_AT_ARP/breadth_coverage/coassembly_AT_to_aer.sam \
  out=~/megahit/megahit_coassembly/aeration_AT_ARP/breadth_coverage/coassembly_AT_to_aer.sam.map.txt

# ARP
# Map reads
bowtie2 -x ~/megahit/megahit_coassembly/aeration_AT_ARP/breadth_coverage/coassembly_aer_contig.index \
  -1 ~/clean_read/ARP*_1.fastq.gz -2 ~/clean_read/ARP*_2.fastq.gz \
  -S ~/megahit/megahit_coassembly/aeration_AT_ARP/breadth_coverage/coassembly_ARP_to_aer.sam.map  -p 16

# Print mapping read number and generate mapping file
$ pileup.sh in=~/megahit/megahit_coassembly/aeration_AT_ARP/breadth_coverage/coassembly_ARP_to_aer.sam \
  out=~/megahit/megahit_coassembly/aeration_AT_ARP/breadth_coverage/coassembly_ARP_to_aer.sam.map.txt
```
## Binning
### [metaWRAP](https://github.com/bxlab/metaWRAP)
#### Installation & environment creation
```
## Install mamba
$ conda install -y mamba
## Download or clone this ripository
$ git clone https://github.com/bxlab/metaWRAP.git

## Configure
$ mkdir ~/db/checkm_db
$ cd ~/db/checkm_db
$ wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
$ tar -xvf *.tar.gz
$ rm *.gz
$ cd
$ checkm data setRoot # Tell CheckM where to find this data before running anything
$ checkm data setRoot db/checkm_db # Tell CheckM where to find this data

## Make metaWRAP executable
$ vi ~/.bash_profile
$ PATH=~/metaWRAP/bin/:$PATH # Add to last row

## Make a new conda environment
$ conda create -y -n metawrap-env python=2.7
$ conda activate metawrap-env

## Install all metaWRAP dependancies
$ conda config --add channels defaults
$ conda config --add channels conda-forge
$ conda config --add channels bioconda
$ conda config --add channels ursky
$ mamba install --only-deps -c ursky metawrap-mg
```
#### metaWRAP binning
```
$ mkdir ~/metawrap_run
$ gunzip ~/clean_read/A*.fastq.gz

# Bin with three different algorithms
$ metawrap binning -o ~/metawrap_run/initial_binning -t 16 -a ~/megahit/megahit_coassembly/aeration_AT_ARP/final.contigs.fa \
  --metabat2 --maxbin2 --concoct ~/clean_read/A*.fastq

# Bin refinement
$ metawrap bin_refinement -o ~/metawrap_run/bin_refinement -t 16 -A ~/metawrap_run/initial_binning/metabat2_bins \
  -B ~/metawrap_run/initial_binning/maxbin2_bins \
  -C  ~/metawrap_run/initial_binning/concoct_bins -c 50 -x 10 -m 56

 # Find the abundaces of the draft genomes (bins) across the samples
 $ metawrap quant_bins -b ~/metawrap_run/bin_refinement/metawrap_50_10_bins -o ~/metawrap_run/bin_quant \
   -a ~/megahit/megahit_coassembly/aeration_AT_ARP/final.contigs.fa ~/clean_read/A*fastq -t 16
```
## Taxonomy Classification of Bin
### [GTDBTk](https://github.com/Ecogenomics/GTDBTk)
Download and alias the GTDB-Tk reference data
[Release database](https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data)
```
# Create the GTDB-Tk environment
$ conda create -n gtdbtk-2.3.0 -c conda-forge -c bioconda gtdbtk=2.3.0
$ conda activate gtdbtk-2.3.0

# Setting environment-specific variables
$ conda env config vars set GTDBTK_DATA_PATH=~/db/gtdbtk/release214

# Check install
$ gtdbtk check_install
$ gtdbtk test --out_dir gtdbtk_test

# Run
$ gtdbtk classify_wf --genome_dir ~/metawrap_run/bin_refinement/metawrap_50_10_bins -x fa \
  --out_dir ~/bins_gtdbdk --scratch_dir scratch.tempory --cpus 16  --skip_ani_screen
```
## Gene Alignment to Bin
### [Prodigal](https://github.com/hyattpd/Prodigal) & [CD-HIT](https://github.com/weizhongli/cdhit)
Gene prediction & reductant sequence remove
```
# Enter environment
$ conda activate prodigal

# ORF perdiction
$ ~/shell_script/prodigal_bin.sh

# Remove reductant sequence
$ ~/shell_script/cdhit_bin.sh
```
### [Diamond](https://github.com/bbuchfink/diamond)
BLAST sequence
```
# Enter environment
$ conda activate diamond

# Run
$ ~/shell_script/diamond_bin.sh
```
### R script after Gene Alignment to bins
```
$ cd ./binning

# QC & merge blast bins with database
Remove 0 bytes dmnd file in ../../airborne_arg_uwtp_result/bins_diamind
$ Rscript bins_ARG_diamond.R
$ Rscript bins_MGE_diamond.R
$ Rscript bins_VF_diamond.R
$ Rscript bin_quality.R
$ Rscript bin_abundance.R
```
### [PROKKA](https://github.com/tseemann/prokka)
```
# Create the prokka environment
$ conda create --name prokka
$ conda install -c conda-forge -c bioconda -c defaults prokka
$ conda activate prokka

# Check install
$ prokka

# Run
$ ~/shell_script/prokka_bin.sh
```