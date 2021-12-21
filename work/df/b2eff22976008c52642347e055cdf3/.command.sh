#!/bin/bash -ue
cooler csort -c1 2 -p1 3 -c2 4 -p2 5 -o test1.pairs.gz test1.valid.pairs ~/22.genome/GRCh38/GRCh38.primary_assembly.genome_no_scaff.fa.sizes
