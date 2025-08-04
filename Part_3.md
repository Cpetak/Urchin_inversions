# Part 3

## Inversion markers

Code for finding inversion markers is described in this notebook:

[finding_markers.ipynb](https://github.com/Cpetak/Urchin_inversions/blob/main/finding_markers.ipynb)

Csv files with SNP inversion markers can be found in the intermediary_files directory.

## GO Enrichment for a specific region

### Step 1 - List of GO terms
Uniprot has GO terms associated to each gene in the urchin genome, here is the link to retrieve that info: https://www.ebi.ac.uk/QuickGO/annotations?taxonId=7668&taxonUsage=exact. 

Used this code to transform uniprot - GO output file into the mapping file topGO expects:

```bash
awk -F "\t" '{print $2"\t"$5}' QuickGO-annotations-1642716310981-20220120.tsv > temp_mapping
sed '$!N; /^\(.*\)\n\1$/!P; D' temp_mapping > temp2_mapping # It deletes duplicate, consecutive lines from a file
awk 'BEGIN{FS="\t"} {for(i=2; i<=NF; i++) { if (!a[$1]) a[$1]=$1FS$i ;else a[$1]=a[$1]","$i};if ($1 != old) b[j++] = a[old];old=$1 } END{for (i=0; i<j; i++) print b[i] }' temp2_mapping > GO_mapping_topGO #it collapses repeated lines into 1, comma separated, output file is in this git repo
```

### Step 2 - List of LOC to Uniprot IDs
Downloaded the gff file from the NCBI genome assembly Spur_5.0 https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002235.5/

Processed file with:

`grep -o 'gene=LOC[0-9]\+' genomic.gff | sed 's/;.*//' | awk -F= '{print $2}' > all_locs_gff`

`sed '$!N; /^\(.*\)\n\1$/!P; D' all_locs_gff > all_locs_gff2`

To get the list of genes (LOC). 32,087 genes.

Uploaded these LOC gene names to https://www.uniprot.org/uploadlists/ From Ensemble Genomes To uniprot. Downloaded results in tsv format. 3 columns: LOC id, uniprot id, uniprot id \_STRPU. Since a LOC id can map to more than 1 uniprot ids, I used the following code to select a specific uniprot id for each LOC to avoid biasing the GO analysis: 

Run: `python loc2uniprot.py`, [loc2uniprot.py](https://github.com/Cpetak/Urchin_inversions/blob/main/loc2uniprot.py)

Output: all_locs_to_uniprotIDs.txt Note: to use, needed to sort file using the `sort all_locs_to_uniprotIDs.txt` command.

### Step 3 - Get list of interesting LOCs and Uniprot IDs
To get the list of genes that are of interest from the NCBI annotation file (Run bash GO_from_region.sh NW_022145594.1 12702886 13424367):

```bash
chr=$1
start=$2
stop=$3

#Most inclusive, includes genes that start before the end of the region (stop) and ends after the start of the region (start). So, it includes genes that span the entire region, the genes that start before the start of the region and end somewhere in the middle, genes that start somewhere in the middle and end outside of the region, and genes that are entirely included.
awk -v mychr=$chr '$1 == mychr' ~/WGS/Urchin_inversions/supp_files/genomic.gff | awk -v mystart=$start '$5 >= mystart {print}' | awk -v mystop=$stop '$4 <= mystop {print}' | awk -v myname='gene' '$3 == myname {print}' | grep -o 'gene=LOC[0-9]\+' | awk -F= '{print $2}' | sort | uniq > ${1}_${2}_${3}_locs1.txt

#More conservative, subset of the first locs list. Only includes genes that start and end within the region.
awk -v mychr=$chr '$1 == mychr' ~/WGS/Urchin_inversions/supp_files/genomic.gff | awk -v mystart=$start '$4 >= mystart {print}' | awk -v mystop=$stop '$5 <= mystop {print}' | awk -v myname='gene' '$3 == myname {print}' | grep -o 'gene=LOC[0-9]\+' | awk -F= '{print $2}' | sort | uniq > ${1}_${2}_${3}_locs2.txt

#Only includes genes that start before the start and end after the start. Includes genes that span entire region. It is a subset of the first locs list. The point is that it lists genes the left breakpoint disrupts.
awk -v mychr=$chr '$1 == mychr' ~/WGS/Urchin_inversions/supp_files/genomic.gff | awk -v mystart=$start '$5 >= mystart {print}' | awk -v mystart=$start '$4 <= mystart {print}' | awk -v myname='gene' '$3 == myname {print}' | grep -o 'gene=LOC[0-9]\+' | awk -F= '{print $2}' | sort | uniq > ${1}_${2}_${3}_locs3_1.txt

#Only includes genes that start before the end and end after the end. Includes genes that span entire region. It is a subset of the first locs list. The point is that it lists genes the right breakpoint disrupts.
awk -v mychr=$chr '$1 == mychr' ~/WGS/Urchin_inversions/supp_files/genomic.gff | awk -v mystop=$stop '$5 >= mystop {print}' | awk -v mystop=$stop '$4 <= mystop {print}' | awk -v myname='gene' '$3 == myname {print}' | grep -o 'gene=LOC[0-9]\+' | awk -F= '{print $2}' | sort | uniq > ${1}_${2}_${3}_locs3_2.txt

#Subset of first locs list, doesn't include genes that are inside the region, not crossing the breakpoints (second locs list). So this set is A - B, where A is all the genes (first locs list), B is only genes in the middle of the region (second locs list)
cat ${1}_${2}_${3}_locs3_1.txt ${1}_${2}_${3}_locs3_2.txt | sort | uniq > ${1}_${2}_${3}_locs3.txt
```
Output: list of gene names (LOC)

Then, to turn LOCS into Uniprot IDs:
```bash
join -t, -1 1 -2 1 -o 1.1,2.2 ${1}_${2}_${3}_locs1.txt ~/WGS/Urchin_inversions/supp_files/all_locs_to_uniprotIDs.txt | awk -F, '{print $2}' > ${1}_${2}_${3}_uniprot1.txt
join -t, -1 1 -2 1 -o 1.1,2.2 ${1}_${2}_${3}_locs2.txt ~/WGS/Urchin_inversions/supp_files/all_locs_to_uniprotIDs.txt | awk -F, '{print $2}' > ${1}_${2}_${3}_uniprot2.txt
join -t, -1 1 -2 1 -o 1.1,2.2 ${1}_${2}_${3}_locs3.txt ~/WGS/Urchin_inversions/supp_files/all_locs_to_uniprotIDs.txt | awk -F, '{print $2}' > ${1}_${2}_${3}_uniprot3.txt
```

### Step 4 - Run GO
Run: `Rscript go_enrichment.R NW_022145594.1 12702886 13424367` (has to be the same inputs as for above steps), [go_enrichment.R](https://github.com/Cpetak/local_pca_pipe/blob/main/go_enrichment.R)

Does GO enrichment for all three kinds of categories listed above for all three kids of GO (biological process, cellular component, molecular function) which are combined into one file. So 3 output files in total.

**All resulting files can be found at: go_enrich_results on Zenodo TODO.**

Code to make the Figure 5 go enrichment subfigure: [go_enrich_figure.R](https://github.com/Cpetak/Urchin_inversions/blob/main/go_enrich_figure.R)

Figure can be found in the intermediary_files directory.

## SnpEff - predicting effect of SNPs

Package description of SnpEff can be found [here](https://pcingola.github.io/SnpEff/snpeff/introduction/).

Downloaded package and tested on example: `java -Xmx8g -jar snpEff.jar -v GRCh37.75 examples/test.chr22.vcf > test.chr22.ann.vcf`

Worked as intended.

Run:

`java -Xmx8g -jar snpEff.jar -v Strongylocentrotus_purpuratus examples/NW_022145594.1_filtered.vcf > 594_test.vcf`

This successfully downloaded the Strongylocentrotus_purpuratus database SnpEff has had already. Knew the name from running: `java -jar snpEff.jar databases`.
Got an error because the chromosome names in the SnpEff database were different from mine so NW_022145594.1 was not found.
In the output text, it listed all its chromosome names and their corresponding lengths, so I just looked up the length of NW_022145594.1 from the NCBI assembly page and search for that in the list. AAGJ06000001.1 has exactly the same number of base pairs, 53101916.

597 - 34141700 - AAGJ06000012.1

600 - 37282239 - AAGJ06000015.1

601 - 35007347 - AAGJ06000016.1

603 - 34285068 - AAGJ06000018.1

606 - 52437917 - AAGJ06000020.1

609 - 39838600 - AAGJ06000003.1

610 - 35917773 - AAGJ06000004.1

So replaced the chromosome name and run:

```bash
sed 's/NW_022145594.1/AAGJ06000001.1/g' examples/NW_022145594.1_filtered.vcf > examples/NW_022145594.1_filtered_renamed.vcf
java -Xmx8g -jar snpEff.jar -v Strongylocentrotus_purpuratus examples/NW_022145594.1_filtered_renamed.vcf > 594_test.vcf
```

Output:
- snpEff_summary.html
- snpEff_genes.txt - table of genes and the variant effects
- 594_test.vcf - additional column in the vcf file, shows the effect of the mutation

Extract data from results files generated above:
cname=594
(all `java -Xmx8g -jar snpEff.jar -v Strongylocentrotus_purpuratus ${chr}_filtered_renamed.vcf > ${cname}_snpeff.vcf`,
`mv snpEff_summary.html ${cname}_snpEff_summary.html`,
`bash get_html_vals.sh ${cname}_snpEff_summary.html ${cname}.csv`)

get_html_vals.sh (that calls extract_intergenic.py) extracts values of interest from the html file.

Repeat with vcf file containing only the inversion:
`bash subset_vcf.sh $chr $actustart $actustop`
cname=594_inv
`java -Xmx8g -jar snpEff.jar -v Strongylocentrotus_purpuratus ${chr}_${actustart}_${actustop}.vcf > ${cname}_snpeff.vcf`,
`mv snpEff_summary.html ${cname}_snpEff_summary.html`,
`bash get_html_vals.sh ${cname}_snpEff_summary.html ${cname}.csv`

**Code above is gathered in get_snpeff_results.sh.**

To get allele frequencies:
cname=594

`bcftools view -i 'ANN[*] ~ "HIGH"' ${cname}_snpeff.vcf > high_${cname}_snpeff.vcf`
`vcftools --vcf high_${cname}_snpeff.vcf --freq --out high_${cname}_snpeff_freq.txt`

`bcftools view -i 'ANN[*] ~ "LOW"' ${cname}_snpeff.vcf > low_${cname}_snpeff.vcf`
`vcftools --vcf low_${cname}_snpeff.vcf --freq --out low_${cname}_snpeff_freq.txt`

`bcftools view -i 'ANN[*] ~ "MODERATE"' ${cname}_snpeff.vcf > mod_${cname}_snpeff.vcf`
`vcftools --vcf mod_${cname}_snpeff.vcf --freq --out mod_${cname}_snpeff_freq.txt`

repeat with cname=594_inv

**Code above is gathered in get_afs.sh**

**Resulting files can be found here: snpeff_results directory.** TODO add to zenodo

**Outputs analysed and figures made in [snpeff.ipynb](https://github.com/Cpetak/Urchin_inversions/blob/main/snpeff.ipynb).**

# Age estimation using GEVA

## Step 1: Set up

Cloned GEVA from: https://github.com/pkalbers/geva?tab=readme-ov-file, then compiled using make.

Got code to run it from: https://github.com/Jcbnunez/Cville-Seasonality-2016-2019/blob/main/CODE/9.0.GEVA_allele_age/2.run_geva_dmel.sh

For phasing, downloaded static release for shapit5 from: https://github.com/odelaneau/shapeit5/releases
Then: `chmod +x phase_common_static`

## Step 2: Making the Guide files

Code to make the guide file:

```python
import sys

def generate_ranges(start, label, stop, filename="guide_file.txt", step=100000):
    start = int(start)  # Ensure start is an integer
    stop = int(stop)  # Ensure stop is an integer

    with open(filename, "w") as f:
        while start + step <= stop:
            end = start + step
            range_col = f"{label}:{start}-{end-1}"
            f.write(f"{start}\t{end}\t{label}\t{step}\t{end-1}\t{start}\t{range_col}\n")
            start = end  # Update start for the next row

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <start> <label> <stop>")
    else:
        start_value = sys.argv[1]
        label_value = sys.argv[2]
        stop_value = sys.argv[3]
        generate_ranges(start_value, label_value, stop_value)
```

To run it: `python make_guide.py 1 NW_022145594.1 500000`

## Step 3: Preparing vcf files and environment

```bash
srun -p bluemoon -N1 --mem=20G --pty bash

spack load bcftools@1.10.2
spack load vcftools@0.1.14
conda activate wgs #has tabix

mychr=NW_022145594.1
cp /users/c/p/cpetak/EG2023/structural_variation/filtered_bcf_files/${mychr}/${mychr}_filtered.vcf .
input_vcf=${mychr}_filtered.vcf

bgzip $input_vcf
tabix ${input_vcf}.gz

./phase_common_static --input ${input_vcf}.gz --region $mychr --output ${input_vcf}_phased.bcf
bcftools view -Ov -o ${input_vcf}_phased.vcf ${input_vcf}_phased.bcf

input_vcf=${input_vcf}_phased.vcf

bgzip $input_vcf
tabix ${input_vcf}.gz

vcftools --gzvcf ${input_vcf}.gz --recode --recode-INFO-all --out ${input_vcf}_processed

bgzip ${input_vcf}_processed.recode.vcf
tabix ${input_vcf}_processed.recode.vcf.gz
```

Put all of the above in a script called prep_files_for_chr.sh, which just takes the chromosome name as an argument and does all of the above. 

Will need to run spack and conda before submitting job!


## Step 4: Running script

Making guide_file: 

```bash
mychr=NW_022145615.1 #guidefile made, not job submitted
input_vcf=${mychr}_filtered.vcf_phased.vcf_processed.recode.vcf.gz
read start stop < <(bcftools view -H ${input_vcf} | awk 'NR==1 {first=$2} {last=$2} END {print first, last}')
python make_guide.py $start $mychr $stop

sbatch --array=1-$( cat guide_file_${mychr}.txt | wc -l) run_geva.sh guide_file_${mychr}.txt 0.01 geva_results $input_vcf
```

To launch jobs one after another use the arraylauncher.py script.

**Analysis and code for the figure can be found in [tmrca.ipynb](https://github.com/Cpetak/Urchin_inversions/blob/main/tmrca.ipynb).**

**Output files can be found in the tmrca_results directory on Zenodo** TODO
