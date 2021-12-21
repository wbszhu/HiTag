#!/bin/bash -ue
echo -n "raw_reads	" >> test1.count
cat test1.raw_reads.count >> test1.count


echo -n "trim_adapter	" >> test1.count
cat test1.trim_adapter.count >> test1.count


echo -n "trim_linker	" >> test1.count
cat test1.trim_linker.count >> test1.count


echo -n "bam2pairs	" >> test1.count
cat test1.bam2pairs.count >> test1.count


echo -n "dedup	" >> test1.count
cat test1.dedup.count >> test1.count


echo -n "UU	" >> test1.count
cat test1.UU.count >> test1.count


echo -n "remove_noise	" >> test1.count
cat test1.remove_noise.count >> test1.count
