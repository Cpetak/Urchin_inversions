#!/bin/bash

# initialize conda
source ~/miniconda3/etc/profile.d/conda.sh

############################################
# PARAMETERS
############################################

filename=NW_022145610.1_30779143_31460853
inv_name=inv9
chr_name=NW_022145610.1
start=30779143
stop=31460853
windowsize=5000

# multiple group types
group_types=("homop" "homoq" "homos" "hetero" "all")

############################################
# LOOP
############################################

for grouptype in "${group_types[@]}"
do

    echo "Running grouptype: $grouptype"

    bash auto_LD_step1.sh $filename $inv_name

    bash auto_LD_step2.sh \
        $chr_name \
        $start \
        $stop \
        $grouptype \
        $inv_name \
        $windowsize \
        FALSE

    conda activate grn

    bash auto_LD_step2.sh \
        $chr_name \
        $start \
        $stop \
        $grouptype \
        $inv_name \
        $windowsize \
        TRUE

    echo "Waiting for SLURM jobs to finish..."

    while [ $(squeue -u cpetak | wc -l) -gt 1 ]
    do
        sleep 30
    done

    echo "Jobs finished."

    Rscript gather_results.R \
        $chr_name \
        $grouptype \
        $windowsize

    Rscript plot_LD.R \
        $chr_name \
        $grouptype \
        $windowsize

    conda deactivate

    rm intermediary_files/${inv_name}/SbatchJob_*

done
