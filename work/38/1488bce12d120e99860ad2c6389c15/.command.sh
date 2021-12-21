#!/bin/bash -ue
pairtools parse -c ~/22.genome/GRCh38/GRCh38.primary_assembly.genome_no_scaff.fa.sizes -o test1.pairs.gz --min-mapq 30 --drop-sam test1.bam
