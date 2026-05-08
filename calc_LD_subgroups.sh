#!/bin/bash

#1 is chromosome name (NW_022145594.1)
#2 is chromosome region (12702886_16793794)
#3 is subgroup (homop)
#4 is LD region calculation start (8702886)
#5 is LD region calculation end (20793794)
#6 is inversion name (inv1)
#7 is num_windows to try (30000)

#Getting subgroup
cd ~/WGS/Urchin_inversions/intermediary_files/$6
filtered=/netfiles/pespenilab_share/urchin_bcfs/$1
#filename=${1}_filtered.vcf_phased.vcf
#filtered=~/WGS/Urchin_inversions/intermediary_files/${6}/${filename}
bcftools view -S ${1}_${2}_vcf_list_${3} -o ${1}_${2}_${3}.vcf ${filtered}/${1}_filtered.vcf
#bcftools view -S ${1}_${2}_vcf_list_${3} -o ${1}_${2}_${3}.vcf ${filename}


#Get specific region
grep '^#' ${1}_${2}_${3}.vcf > ${1}_${2}_${3}_header
mystart=$4
myend=$5
grep -v \# ${1}_${2}_${3}.vcf | awk -v myvariable=$mystart '$2 >= myvariable' | awk -v myvariable=$myend '$2 <= myvariable' > temp.vcf
outfilename=${1}_${4}_${5}_${3}.vcf
cat ${1}_${2}_${3}_header temp.vcf > $outfilename
rm temp.vcf

