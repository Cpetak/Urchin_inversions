# Input data processing

The pipeline will start from filtered bcf files, one for each of the 21
chromosomes.

Below are the steps that I followed to get from the raw sequence files to
these bcf files.

## Mapped reads to the reference

### The reference genome

<https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002235.5/>

N50: 37.3 Mb

### The mapping algorithm

Input: List of read files (R1 and R2)

```bash
while read line ; do
        F1=$(cut -d ' ' -f1 <<< $line)
        F2=$(cut -d ' ' -f2 <<< $line)
        echo "$F1 -- $F2"
        FILE=$(mktemp)
        cat header.txt >> $FILE
        echo "spack load samtools@1.10" >> $FILE
        echo "spack load bwa@0.7.17" >> $FILE
        ref="/users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna"
        out_name=$(cut -d '.' -f1 <<< $F1)
        echo "bwa mem -t 1 -M $ref /users/c/p/cpetak/WGS/all_fastqs/$F1 /users/c/p/cpetak/WGS/all_fastqs/$F2 | samtools view -S -b > /users/c/p/cpetak/WGS/BWA_out/$out_name.bam" >> $FILE
        sbatch $FILE
        sleep 0.5
        rm $FILE
done < $1
```

The Burrows-Wheeler Alignment Tool (BWA) MEM algorithm was used for
mapping the raw reads to the S. purpuratus reference genome (Spur ver.
5.0, scaffold N50 ∼37 Mbp). The average coverage for each individual was
6.42±0.78, with an average mapping rate of 81.6±0.01.

## Called variants for each chromosome across all individuals

Input:

-   21 chromosome names
-   list_of_files.txt, 140 lines, line 1:
    `/users/c/p/cpetak/WGS/BWA_out/BOD_18170X61_200925_A00421_0244_AHKML5DSXY_S81_L002_R1_001.rmdup.bam`

```bash
while read line ; do
	echo "$line"
	FILE=$(mktemp)
  cat header.txt >> $FILE
  ref="/users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna"
  echo "echo "${line}" " >> $FILE
  echo "bcftools mpileup -r $line -f $ref --bam-list list_of_files.txt | bcftools call -mv -Ob -o multi_bam_${line}.bcf" >> $FILE
  sbatch $FILE
  sleep 0.5
  rm $FILE
done < $1
```

## Filtering the bcf files

```bash
#!/bin/sh

mychr="NW_022145594.1"
myfolder="/users/c/p/cpetak/EG2023/structural_variation/bcf_files"
suppfolder="/users/c/p/cpetak/WGS/local_pca_pipe/supp_files"

#pre filter: 10.854.890

bcftools view -e 'QUAL <= 40 || DP < 560 || MQB < -3 || RPB < -3 || RPB > 3 || AN < 238' ${myfolder}/multi_bam_${mychr}.bcf > ${myfolder}/${mychr}_filtered.vcf
#note: after each step, the output is vcf, which needs to be converted into bcf
bcftools view -Ob ${myfolder}/${mychr}_filtered.vcf > ${myfolder}/${mychr}_filtered.bcf

#post first filter (aka old filtering): 3.999.255

# taking out 3 outliers: CAP_18170X101, FOG_18170X127 and FOG_18170X128
bcftools view -S ${suppfolder}/all_rmdups_noout.txt -o ${myfolder}/${mychr}_filtered_noout.vcf ${myfolder}/${mychr}_filtered.bcf
# Allele count in genotypes for each ALT allele, keep only SNPs, and only biallelic:
bcftools view -e 'AC < 14' --exclude-types 'indels,mnps,ref,bnd,other' -m2 -M2 ${myfolder}/${mychr}_filtered.bcf > ${myfolder}/${mychr}_filtered_noout.vcf
# remove repetitive elements:
#found a paper that used RepeatModeler 2.0.1 to find repetitive regions in the Spur assembly (same as what we are using here) - took their dataset and removed any variations that fell within any of those regions. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8615465/#:~:text=Repetitive%20elements%20(REs)%20occupy%20a,of%20biological%20processes%20remains%20unknown

bcftools view -T ^${suppfolder}/Spur_repeats_02 ${myfolder}/${mychr}_filtered.bcf -o ${myfolder}/${mychr}_filtered_noout.vcf

#post second filter: 402.500

# Convert the filtered vcf into the bcf file type which is the type the R package will be expecting

bcftools view -Ob ${myfolder}/${mychr}_filtered.vcf > ${myfolder}/${mychr}_filtered.bcf

# Index the filtered bcf file. This will make the file more searchable by the algorythm reading it.

bcftools index ${myfolder}/${mychr}_filtered.bcf
```

Output: bcf files used to be in the
`/users/c/p/cpetak/EG2023/structural_variation/filtered_bcf_files`
directory.

Now they are in `/netfiles/pespenilab_share/urchin_bcfs`.

# Running lostruct

On the github page, they provide an Rscript to show how to use the package. I copied this:
[run_lostruct.R](https://github.com/Cpetak/Urchin_inversions/blob/main/run_lostruct.R)

No changes, except output folder has chromosome name instead of random id.

```bash     
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
```

To run: `bash 2_1_local_pca.sh NW_022145594.1`

## Visualise and gather data

They also have an example of plotting the data, I copied this as well.
[summarize_run.Rmd](https://github.com/Cpetak/Urchin_inversions/blob/main/summarize_run.Rmd)

Changes:

-   do.pdfs \<- TRUE in line 14,

-   commented out warning in line 78, breaks code otherwise

-   added saving of corner pcas

-   added getting the percent corner (i.e. alpha) as a variable read
    from a file

Once the 2.1 script finishes, you should have a folder called lostruct_results and a folder in it, named type_snp_size_10000_chromosome_NW_022145594.1.

```bash
#conda activate grn
#spack load bcftools@1.10.2

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
```

To run: `bash 2_2_local_pca_vis.sh NW_022145594.1 snp 10000 0.01` It is chromosome, type, size, alpha (percent outliers considered). 
To run for every snp length independently, a one-line:
`for snp in 500 1000 5000 10000; do bash 2_2_local_pca_vis.sh NW_022145594.1 snp $snp 0.05; done`

Then can rerun specific ones to adjust alpha.

**Results of the local PCA analysis can be found here: TODO Zenodo link. Already uploaded, not published yet.**

## Finding putative breakpoints

Automated way of finding the genomic coordinates of the outlier regions.

INPUT: MDS output files for a specific chromosome, for a specific type of window size (snp vs bp)

Python file, conda activate wgs, [finding_outlier_windows.py](https://github.com/Cpetak/Urchin_inversions/blob/main/finding_outlier_windows.py)

To run: `python finding_outlier_windows.py NW_022145594.1 snp`

You can add an optional argument `--thr 0.2` to adjust the MDS outlier threshold.

OUTPUT: Outlier region genomic coordinates for a specific chromosome for each window size tested.

E.g., for NW_022145594.1, two spikes, (from 500 snp window data, all filtering): 

12702886-16793794 whole region, first spike: 12702886 - 13424367, second spike: 15422748-16793794

**Resulting breakpoint coordinates are summarized here: [breakpoints.txt](https://github.com/Cpetak/Urchin_inversions/blob/main/breakpoints.txt)**


