#!/bin/bash -ue
python /public/home/luzhang/00.software/HiTag/./scripts/format_fithichip_config.py test1.validPairs /public/home/luzhang/00.software/HiTag/./test_data/H3K27ac_K562_pup5_qup2_peak.bed test1.fithichip.out 5000 ~/22.genome/GRCh38/GRCh38.primary_assembly.genome_no_scaff.fa.sizes 0 > test1_fithichip.config
FitHiChIP_HiCPro.sh -C test1_fithichip.config
