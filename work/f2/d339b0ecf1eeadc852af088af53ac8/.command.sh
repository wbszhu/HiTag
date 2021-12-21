#!/bin/bash -ue
echo "--------------------------" >> test2.stat
echo dedup_stat >> test2.stat
echo "--------------------------" >> test2.stat
cat test2.dedup.stats.txt >> test2.stat
echo >> test2.stat


echo "--------------------------" >> test2.stat
echo UU_stats >> test2.stat
echo "--------------------------" >> test2.stat
cat test2.UU.stats.txt >> test2.stat
echo >> test2.stat


echo "--------------------------" >> test2.stat
echo remove_noise >> test2.stat
echo "--------------------------" >> test2.stat
cat test2.noise_reduce.log >> test2.stat
echo >> test2.stat
