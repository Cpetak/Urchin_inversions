#!/usr/bin/env bash
#SBATCH -J BayPassrun 
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 30:00:00 
#SBATCH --mem 60G 
#SBATCH -p general
#SBATCH --mail-type=END
#SBATCH --mail-user=daniel.sadler@uvm.edu
#SBATCH -e errors_%A_%a.txt
#SBATCH -o output_%A_%a.txt
#SBATCH --array=1-1000 

cd /users/d/s/dsadler1/scratch/SeaUrchinEnvPop/GDS/BayPass/random_chunks

module load gcc/13.3.0
module load openmpi/5.0.5

# Split genotype file into 1000 chunks

# Run BayPass for each chunk
/users/d/s/dsadler1/programmes/baypass_public-master/sources/./g_baypass -gfile genotypefilechunk_$SLURM_ARRAY_TASK_ID.txt -outprefix full_$SLURM_ARRAY_TASK_ID -nthreads 8