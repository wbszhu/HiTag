#!/bin/bash -ue
l=$(zcat test1_R1_val_1.fq.gz | wc -l | awk '{print $1}')
echo $(($l / 4)) > test1.trim_adapter.count
