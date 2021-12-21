#!/bin/bash -ue
trim_galore -q 20 --phred33 --paired --illumina --trim-n --gzip test1_R1.fq.gz test1_R2.fq.gz
