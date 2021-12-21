#!/bin/bash -ue
l=$(zcat test2_R1_val_1.fq.gz | wc -l | awk '{print $1}')
echo $(($l / 4)) > test2.trim_adapter.count
