# Coverage

### Step 1

A) Make a file called myposi.bed with one line in it: `NW_022145594.1 15422748 16793794`, chromosome, start, stop of region to get coverage for.

B) Make a file called bam_list: `/netfiles/pespenilab_share/urchin_bams/BOD_18170X61_200925_A00421_0244_AHKML5DSXY_S81_L002_R1_001.rmdup.bam`, etc, listing all bam files.

C) Make 3 files: hetero_list_cov, homop_list_cov, homoq_list_cov, listing names of coverage files for each group.

### Step 2

Calculate coverage for each individual, for that specific region.

```bash
while read line ; do
        echo "$line"
        FILE=$(mktemp)
        pos_file='~/WGS/local_pca_pipe/myposi.bed'
        cat header.txt >> $FILE
        echo "spack load samtools@1.10" >> $FILE
        out_name=$(cut -d '.' -f1 <<< $line)
        out_name2=$(echo "$out_name" | cut -d'/' -f5)
        echo "samtools depth -b $pos_file $line > ~/WGS/inversion_results/${out_name2}.coverage" >> $FILE
        sbatch $FILE
        #cat $FILE
        sleep 0.1
        rm $FILE
done < $1
```
Run: `bash cov_launcher.sh bam_list`

### Step 3

Move .coverage files in the 2 folders: left_bp (breakpoint area), right_bp.

Run: `python coverage_plotting.py`, [coverage_plotting.py](https://github.com/Cpetak/Urchin_inversions/blob/main/coverage_plotting.py)

