#!/bin/bash -ue
bwa mem -SP5M -t 20 ~/22.genome/GRCh38/bwa_index/GRCh38.primary_assembly.genome_no_scaff.fa test1_1.valid.fastq test1_2.valid.fastq | samtools view -b > test1.bam
