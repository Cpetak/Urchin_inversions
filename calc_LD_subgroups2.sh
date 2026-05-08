#!/bin/bash

#1 is chromosome name (NW_022145594.1)
#2 is chromosome region (12702886_16793794)
#3 is subgroup (homop)
#4 is LD region calculation start (8702886)
#5 is LD region calculation end (20793794)
#6 is inversion name (inv1)
#7 is num_windows to try (30000)

cd ~/WGS/Urchin_inversions/intermediary_files/$6

Rscript ../../vcf2gds.R ${1}_${4}_${5}_${3}.vcf ${1}_${4}_${5}_${3}
Rscript ../../LD_get_num_windows.R ~/WGS/Urchin_inversions/intermediary_files/$6/${1}_${4}_${5}_${3}.gds $7
