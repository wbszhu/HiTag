#!/bin/bash -ue
echo -n "raw_reads	" >> test2.count
cat test2.raw_reads.count >> test2.count


echo -n "trim_adapter	" >> test2.count
cat test2.trim_adapter.count >> test2.count


echo -n "trim_linker	" >> test2.count
cat test2.trim_linker.count >> test2.count


echo -n "bam2pairs	" >> test2.count
cat test2.bam2pairs.count >> test2.count


echo -n "dedup	" >> test2.count
cat test2.dedup.count >> test2.count


echo -n "UU	" >> test2.count
cat test2.UU.count >> test2.count


echo -n "remove_noise	" >> test2.count
cat test2.remove_noise.count >> test2.count
