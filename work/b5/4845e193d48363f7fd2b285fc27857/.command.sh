#!/bin/bash -ue
python /public/home/luzhang/00.software/HiTag/./scripts/format_fithichip_config.py test2.validPairs /public/home/luzhang/00.software/HiTag/./test_data/H3K27ac_K562_pup5_qup2_peak.bed test2.fithichip.out 5000 ~/22.genome/GRCh38/GRCh38.primary_assembly.genome_no_scaff.fa.sizes 0 > test2_fithichip.config
FitHiChIP_HiCPro.sh -C test2_fithichip.config
