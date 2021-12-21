#!/bin/bash -ue
l=$(cat test1_1.valid.fastq | wc -l | awk '{print $1}')
echo $(($l / 4)) > test1.trim_linker.count
