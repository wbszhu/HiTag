#!/bin/bash -ue
cat test2.valid.pairs | grep -v "^#" | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$6,$4,$5,$7}' > test2.validPairs
