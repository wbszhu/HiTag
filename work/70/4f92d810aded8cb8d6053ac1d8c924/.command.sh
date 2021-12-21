#!/bin/bash -ue
echo "--------------------------" >> test1.stat
echo dedup_stat >> test1.stat
echo "--------------------------" >> test1.stat
cat test1.dedup.stats.txt >> test1.stat
echo >> test1.stat


echo "--------------------------" >> test1.stat
echo UU_stats >> test1.stat
echo "--------------------------" >> test1.stat
cat test1.UU.stats.txt >> test1.stat
echo >> test1.stat


echo "--------------------------" >> test1.stat
echo remove_noise >> test1.stat
echo "--------------------------" >> test1.stat
cat test1.noise_reduce.log >> test1.stat
echo >> test1.stat
