# Airborne ARGs in UWTP
A comprehensive metagenomic pipeline for airborne ARGs (Antibiotic Resistance Genes) in UWTP (Urban wastewater Treatment Plant)

## Table of Contents
1. [Anaconda Installation](https://github.com/yenjh3910/Metagenomic_pipeline/blob/main/README.md#anaconda-installation)
2. [Taxanomic Profile](https://github.com/yenjh3910/Metagenomic_pipeline/blob/main/README.md#taxanomic-profile)
3. [ARGs Profile](https://github.com/yenjh3910/Metagenomic_pipeline/blob/main/README.md#args-profile)
4. [Functional Profile](https://github.com/yenjh3910/Metagenomic_pipeline/blob/main/README.md#functional-profile)

## Anaconda Installation
https://www.anaconda.com/products/distribution
```
$ sudo apt update && upgrade
$ sudo apt install curl
$ mkdir tmp
$ cd tmp

# Install
$ curl -O https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh
$ bash Anaconda3-2022.10-Linux-x86_64.sh
$ source ~/.bashrc
$ conda info

$ cd
$ rm -r tmp
```
### Bioconda Installation
https://bioconda.github.io/
```
# Bioconda set up channels
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
$ conda config --set channel_priority strict
# View channel
$ conda config --show channels
```
## Taxanomic Profile
### Kraken2
https://github.com/DerrickWood/kraken2/wiki/Manual  
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
cd ~/db
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
### Bracken
https://github.com/jenniferlu717/Bracken
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
## ARGs Profile
### ARGs-OAP
https://github.com/xinehc/args_oap
#### Installation
```
$ conda create -n args_oap
$ conda activate args_oap
$ conda config --set channel_priority flexible
$ conda install -c bioconda -c conda-forge args_oap
```
#### Usage
```
# Stage one
$ args_oap stage_one -i ~/clean_read -o ~/args_oap/stage_one_output -f fastq -t 16

# Stage two (e_value: 1e-7, identity: 80, aa_length,25)
$ args_oap stage_two -i ~/args_oap/stage_one_output -o ~/args_oap/stage_two_output -t 16
```
## Functional Profile
### HUMAnN 3.0
https://github.com/biobakery/humann/tree/8d69f3c84ca7bfd7519ced7fcf94b8356c915090
#### Installation
```
$ conda create --name mpa -c bioconda python=3.7
$ conda activate mpa
$ conda install humann
```
