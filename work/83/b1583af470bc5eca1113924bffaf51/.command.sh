#!/bin/bash -ue
pairtools dedup --mark-dups --output-stats test2.dedup.stats.txt -o test2.dedup.pairs.gz test2.sorted.pairs.gz
