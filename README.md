# HiTag nextflow pipeline

Gene expression is not only influenced by linear regulatory elements, but an increasing number of studies have shown that chromatin spatial interactions also play an important role in the regulation of gene expression. The discovery of this phenomenon is based on the continuous advances in chromatin conformation capture technology.

More than a decade has passed since the first chromosomal conformation capture 3C technology became available. During this time, the quality and quantity of capture techniques from chromatin space have increased significantly. Techniques such as circular chromosome conformation capture (4C), 3C-carbon copies derived from 3C (5C) have opened up the possibility of studying long-distance interactions between chromatin. The advent of Hi-C was a major breakthrough in the study of the spatial structure of biological genomes, and Hi-ChIP and ChIA-PET were an important step towards the detection of specific protein-mediated spatial interactions. 

Hi-C, ChIA-PET, HiChIP and PLAC-seq require a largenumber of cells and are cumbersome to use in certain cell-sparse samples due to their low capture efficiency. chromatin conformation, we developed the Hi-Tag.
## Dependencies

The codebase relies on the following dependancies (tested version provided in 
parentheses):

```
trim-galore
pysam
bwa
pairtools
bedtools
cooler
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

### install pysam

```bash
$ conda config --add channels r
$ conda config --add channels bioconda
$ conda install pysam
```
OR you can install through pypi
```
pip install pysam
```
### install fithichip
Find details in this way https://ay-lab.github.io/FitHiChIP/usage/installation.html

## config
When you run the nextflow or fithichip commands and need to generate your own config file quickly and easily, you may want to
config for nextflow
```
python scripts/nextflow_config.py linkerA linkerB min_mapq gsize bwa_index Enzyme.bed> nextflow.config
```
config for fithichip
```
python scripts/format_fithichip_config.py ValidPairs PeakFile OutDir LowDistThr genome_size params.use_P2P_background > fithichip.config
```

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
# run nextflow
$ nextflow run main.nf
```
## Note
Here, linkerA and linkerB represents the sense-strand(F) and antisense-strand(R) of the bridge linker, respectively.
Like this
```
linker A (F)
-------->**ACGCGATATCTTATC**TGACT-------
---------TGCGC**TATCGAATAGACTGA**<------
                              linker B(R)
```
