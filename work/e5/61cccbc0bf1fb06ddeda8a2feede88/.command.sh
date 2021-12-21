#!/bin/bash -ue
trimLinker -e 2 -t 12 -m 1 -k 2 -l 16 -o . -n test2 -A ACGCGATATCTTATC -B AGTCAGATAAGATAT test2_R1_val_1.fq.gz test2_R2_val_2.fq.gz
