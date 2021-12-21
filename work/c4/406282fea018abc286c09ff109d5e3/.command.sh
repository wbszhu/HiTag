#!/bin/bash -ue
l=$(cat test2_1.valid.fastq | wc -l | awk '{print $1}')
echo $(($l / 4)) > test2.trim_linker.count
