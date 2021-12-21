#!/bin/bash -ue
trimLinker -e 2 -t 12 -m 1 -k 2 -l 16 -o . -n test1 -A ACGCGATATCTTATC -B AGTCAGATAAGATAT test1_R1_val_1.fq.gz test1_R2_val_2.fq.gz
