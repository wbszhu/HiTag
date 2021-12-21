#!/bin/bash -ue
trim_galore -q 20 --phred33 --paired --illumina --trim-n --gzip test2_R1.fq.gz test2_R2.fq.gz
