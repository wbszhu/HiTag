#!/bin/bash -ue
pairtools dedup --mark-dups --output-stats test1.dedup.stats.txt -o test1.dedup.pairs.gz test1.sorted.pairs.gz
