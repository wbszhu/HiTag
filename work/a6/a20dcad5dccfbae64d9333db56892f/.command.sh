#!/bin/bash -ue
l=$(zcat test1_R1.fq.gz | wc -l | awk '{print $1}')
echo $(($l / 4)) > test1.raw_reads.count
