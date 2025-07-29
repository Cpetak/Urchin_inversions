#Script to get subsetted bcf file

chrom=$1
mystart=$2
myend=$3

input_vcf="/users/c/p/cpetak/scratch/EG2023/structural_variation/filtered_bcf_files/${chrom}/${chrom}_filtered.vcf"

grep -v \# $input_vcf | awk -v myvariable=$mystart '$2 >= myvariable' | awk -v myvariable=$myend '$2 <= myvariable' > temp.vcf

outfilename=${1}_${2}_${3}.vcf
cat vcf_header_noout temp.vcf > $outfilename

rm temp.vcf
