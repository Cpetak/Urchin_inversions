#!/bin/bash

#NW_022145594.1_12670717_16440127_homop_nums.txt

#$1 NW_022145594.1_12670717_16440127
#$2 num clusters

file_path="vcf_col_ids"

if [ $2 == 3 ]; then

file_homop="${1}_homop_nums.txt"
file_homoq="${1}_homoq_nums.txt"

homop=()
while IFS= read -r line; do
  homop+=("$line")
done < $file_homop

homoq=()
while IFS= read -r line; do
      homoq+=("$line")
    done < $file_homoq

for row_number in "${homop[@]}"; do
  sed -n "${row_number}p" "$file_path" >> "${1}_vcf_list_homop"
done

for row_number in "${homoq[@]}"; do
    sed -n "${row_number}p" "$file_path" >> "${1}_vcf_list_homoq"
  done

else

for c in $(seq 0 $2)
do
  echo "number $c"
  myfile="${1}_clust${c}.txt"
  geno=()
  while IFS= read -r line; do
    geno+=("$line")
  done < $myfile
  for row_number in "${geno[@]}"; do
    sed -n "${row_number}p" "$file_path" >> "${1}_vcf_list_clust${c}"
  done

done

fi
