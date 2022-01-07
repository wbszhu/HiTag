# HiTag nextflow pipeline

Gene expression is not only influenced by linear regulatory elements, but an increasing number of studies have shown that chromatin spatial interactions also play an important role in the regulation of gene expression. The discovery of this phenomenon is based on the continuous advances in chromatin conformation capture technology.

More than a decade has passed since the first chromosomal conformation capture 3C technology became available. During this time, the quality and quantity of capture techniques from chromatin space have increased significantly. Techniques such as circular chromosome conformation capture (4C), 3C-carbon copies derived from 3C (5C) have opened up the possibility of studying long-distance interactions between chromatin. The advent of Hi-C was a major breakthrough in the study of the spatial structure of biological genomes, and Hi-ChIP and ChIA-PET were an important step towards the detection of specific protein-mediated spatial interactions. 

Hi-C, ChIA-PET, HiChIP and PLAC-seq require a largenumber of cells and are cumbersome to use in certain cell-sparse samples due to their low capture efficiency. chromatin conformation, we developed the Hi-Tag.
## Dependencies

The codebase relies on the following dependancies (tested version provided in 
parentheses):

```
trim-galore
bwa
pairtools
bedtools
cooler
ChIA-PET2
Fithichip
HiC-Pro

```

## Description of the HiTag

Compared with the traditional HiChIP and ChIA-PET, the experimental cycle time is greatly shortened and the samples to be sequenced can be obtained in only two days. The present invention can obtain the same chromatin spatial structure map as the traditional technique with relatively small amount of data, which greatly saves the sequencing cost. The data signal-to-noise ratio of this invention is higher than that of the traditional method, which effectively increases the proportion of valid data in the data.


## Installation

conda create a clean environment

```bash
$ conda env create -f HiTag.yml  -p ~/.conda/envs/HiTag
$ conda activate HiTag
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
$ conda activate HiTag
# add paths to $PATH
$ export PATH=/path/to/nextflow:$PATH
$ export PATH=/path/to/HiC-Pro:$PATH
$ export PATH=/path/to/FitHiChIP_HiCPro.sh:$PATH
$ export PATH=/path/to/trimLinker:$PATH
# run nextflow
$ nextflow run main.nf
```
