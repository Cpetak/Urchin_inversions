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

**Results of the local PCA analysis can be found here: TODO Zenodo link.**

Note: Figures and code to make supp figure 1 with the "Genome-wide patterns of local PCA" can be found in the intermediary_files directory and in the global_local_PCA.ipynb.

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

## Grouping individuals by genotype for specific inversion

INPUT: Coordinates of region of interest, bcf file of chromosome

Now that we have regions of interest (to recap: local pca -\> mds values
-\> outlier region coordinates) we should first look at the PCA of the
individuals based on that region. Yes, technically, this is what local
PCA is doing already, but it is cleaner and safer to start fresh from
the bcf files and do our own, independent PCA. The result should be the
same or similar! it might not be exactly the same because the PCAs that
are plotted as part of the local PCA pipe line only include outlier
windows, whereas here we are going to take all SNPs in a specific region
for the PCA analysis, and this might include windows that were not
outliers in the local PCA analysis.

For all of the steps below, conda activate wgs. (grn is for local PCA and LD only)

### Step 1: Make a vcf file that only includes the region of interest for our chromosome.

```bash
chrom=$1
mystart=$2
myend=$3

input_vcf="/users/c/p/cpetak/EG2023/structural_variation/filtered_bcf_files/${chrom}/${chrom}_filtered.vcf"

grep -v \# $input_vcf | awk -v myvariable=$mystart '$2 >= myvariable' | awk -v myvariable=$myend '$2 <= myvariable' > temp.vcf

outfilename=${1}_${2}_${3}.vcf
cat vcf_header_noout temp.vcf > $outfilename

rm temp.vcf
```

run: `bash subset_bcf.sh NW_022145594.1 12670717 16440127` where it is
chromosome, start, end. Looks for a vcf file in the common data folder
for that chromosome, output is also a vcf in the curr directory

### Step 2: Convert vcf to gds

```bash
library(SeqArray)
library(SNPRelate)

seqParallelSetup(cluster=10, verbose=TRUE)

args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
out_file <- args[2]

print(args)
#print(paste(args[1],"sometime"))

snpgdsVCF2GDS(vcf_file,paste(out_file,".gds",sep=""),verbose=T)
```

run:
`Rscript vcf2gds.R NW_022145594.1_12670717_16440127.vcf NW_022145594.1_12670717_16440127`
input filename, output file name (it will add .gds)

### Step 3: Do PCA from the gds file

```R
library("SNPRelate")

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
genofile <- snpgdsOpen(filename)

ccm_pca<-snpgdsPCA(genofile, autosome.only=FALSE)
dim1<-ccm_pca$eigenvect[,1]
dim2<-ccm_pca$eigenvect[,2]

eigenval_list <- ccm_pca$eigenval
eigenval_list[is.nan(eigenval_list)] <- 0
esum <- sum(eigenval_list)
e1 <- eigenval_list[1] / esum * 100
e2 <- eigenval_list[2] / esum * 100

my_df1 <- as.data.frame(dim1)
my_df2 <- as.data.frame(dim2)

special_character <- "\\."
split_string <- strsplit(filename, special_character)[[1]]
result1 <- paste(paste(head(split_string, -1), collapse = "."),"_dim1.csv",sep="")
result2 <- paste(paste(head(split_string, -1), collapse= "."),"_dim2.csv",sep="")

eresult <- paste(paste(head(split_string, -1), collapse= "."),"_perc_explained.csv",sep="")

writeLines(c(as.character(e1), as.character(e2)), eresult)

write.csv(my_df1, result1, row.names=FALSE)
write.csv(my_df2, result2, row.names=FALSE)
```

run: `Rscript do_pca.R NW_022145594.1_12670717_16440127.gds` does PCA
and writes PC1 and PC2 to csv files.

### Step 4: Assign individuals to genotype groups, make plots, correlate with latitude

Python file, [genotype_by_PCA.py](https://github.com/Cpetak/Urchin_inversions/blob/main/genotype_by_PCA.py)

run:
`python genotype_by_PCA.py NW_022145594.1_12702886_16793794_dim1.csv NW_022145594.1_12702886_16793794_dim2.csv 3 NW_022145594.1_12702886_16793794_perc_explained.csv`
3 is number of clusters, output is 2 csv files with individual IDs
(1-140) for each homozygote group and PCA, elbow, pie, line and map
plots, and chi-square results to txt file.

OUT: AA, Aa or aa for each individual and PCA plot, with map

**Resulting genotype frequencies are summarized here: [genotype_freqs.txt](https://github.com/Cpetak/Urchin_inversions/blob/main/genotype_freqs.txt)**

**Resulting plots are in this directory: intermediary_files**

**For changing/recreating figure 1: [Fig1.md](https://github.com/Cpetak/Urchin_inversions/blob/main/Fig1.md)**

**For changing/recreating subplots with PCA where individuals are colored by population: [PCA_notebook.ipynb](https://github.com/Cpetak/Urchin_inversions/blob/main/PCA_notebook.ipynb)**

Note: PCA_notebook.ipynb also has code for 6-way grouping for inversion 5.
