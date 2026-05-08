#!/bin/bash

#1 is inversion label (inv1)
#2 is subgroup type (homop)
#3 is chromosome name (NW_022145594.1)
#4 is LD region calculation start (8702886)
#5 is LD region calculation end (20793794)
#6 is num_windows (30000)

cd ~/WGS/Urchin_inversions
mkdir makegrid_tempfiles_${3}_${2}_${6}
cd ~/WGS/Urchin_inversions/intermediary_files/$1

sbatch ../../makegrid_launcher.sh ~/WGS/Urchin_inversions/intermediary_files/$1/${3}_${4}_${5}_${2}.gds $6
