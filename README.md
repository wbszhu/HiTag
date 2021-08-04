# BL-HiChIP nextflow pipeline

## Installation
conda create a clean enviroment
```bash

$conda create -p ~/.conda/envs/bl-hichip
$conda activate bl-hichip
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


## Run

```bash
# edit config file
$ vim nextflow.config
# enter environment
$ conda activate bl-hichip
# add paths to $PATH
$ export PATH=/path/to/nextflow:$PATH
$ export PATH=/path/to/HiC-Pro:$PATH
$ export PATH=/path/to/FitHiChIP_HiCPro.sh:$PATH
$ export PATH=/path/to/trimLinker:$PATH
# run nextflow
$ nextflow run main.nf
```
