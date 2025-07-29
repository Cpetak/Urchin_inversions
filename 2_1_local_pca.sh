
#conda activate grn
#spack load bcftools@1.10.2

echo "Chromosome number: $1"

input_dir="/users/c/p/cpetak/EG2023/structural_variation/filtered_bcf_files/${1}"

echo $input_dir

for snp in 500 1000 5000 10000
do
    echo "$snp"
    FILE=$(mktemp)
    cat header.txt >> $FILE
    echo "Rscript ~/WGS/local_pca_pipe/run_lostruct.R -i ${input_dir} -t snp -s ${snp} -I ~/WGS/local_pca_pipe/sample_info.tsv -c ${1}" >> $FILE
    sbatch $FILE
    #cat $FILE
    sleep 0.1
    rm $FILE
done
