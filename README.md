# DRP-DNA-DNA nextflow pipeline

## Installation
conda create a clean enviroment
```bash

$ conda create -p ~/.conda/envs/DRP-DNA-DNA
$ conda activate DRP-DNA-DNA
```
### conda 

```bash
$ conda install -c bioconda trim-galore
$ conda install -c bioconda bwa
$ conda install -c bioconda pairtools
$ conda install -c bioconda bedtools 
$ conda install -c bioconda cooler
```
### install trimLinker

```bash

$ git clone https://github.com/GuipengLi/ChIA-PET2.git
$ tar -zxvf ChIA-PET2.tar.gz
$ cd ChIA-PET2
$ chmod +x bin/ChIA-PET2
$ make

```
you can get trimLinker path = ~/00.software/ChIA-PET2/bin/trimLinker

### install fithichip
Find details in this way https://ay-lab.github.io/FitHiChIP/usage/installation.html

## Run

```bash
# edit config file
$ vim nextflow.config
# enter environment
$ conda activate DRP-DNA-DNA
# add paths to $PATH
$ export PATH=/path/to/nextflow:$PATH
$ export PATH=/path/to/HiC-Pro:$PATH
$ export PATH=/path/to/FitHiChIP_HiCPro.sh:$PATH
$ export PATH=/path/to/trimLinker:$PATH
# run nextflow
$ nextflow run main.nf
```
