
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
