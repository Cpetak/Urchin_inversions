
#!/bin/bash

#1 is filename, eg NW_022145594.1_12702886_16793794
#2 is inversion folder

cd ~/WGS/Urchin_inversions/intermediary_files/${2}

grep -Fvx -f ${1}_vcf_list_homoq ../../vcf_col_ids > temp_hetero
grep -Fvx -f ${1}_vcf_list_homop temp_hetero > temp_hetero2
mv temp_hetero2 ${1}_vcf_list_hetero
rm temp_hetero

cat ${1}_vcf_list_homoq ${1}_vcf_list_homop > ${1}_vcf_list_homos

cat ${1}_vcf_list_homos ${1}_vcf_list_hetero > ${1}_vcf_list_all
