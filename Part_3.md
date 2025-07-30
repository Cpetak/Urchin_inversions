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

#ITT TARTOK

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

