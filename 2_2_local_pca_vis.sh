
echo "Chromosome number: $1"
echo "Type: $2"
echo "Size: $3"
echo "Percent corner: $4"

input_dir="/users/c/p/cpetak/EG2023/structural_variation/filtered_bcf_files/${1}"

echo $input_dir

cd ~/WGS/inversion_results/lostruct_results/type_${2}_size_${3}_chromosome_${1}
echo $4 > percent_file.txt #the R script below will be looking for this file! easier than figuring out how to pass in as an argument

cd ~/WGS/inversion_results

FILE=$(mktemp)
cat header.txt >> $FILE
echo "Rscript -e 'templater::render_template(\"~/WGS/local_pca_pipe/summarize_run.Rmd\",output=\"~/WGS/inversion_results/lostruct_results/type_${2}_size_${3}_chromosome_${1}/run_summary.html\",change.rootdir=TRUE)'" >> $FILE
sbatch $FILE
#cat $FILE
sleep 0.1
rm $FILE
