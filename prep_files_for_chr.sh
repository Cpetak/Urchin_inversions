#!/bin/sh

#Specify a partition
#SBATCH --partition=bluemoon
# Request nodes
#SBATCH --nodes=1
# Request some processor cores
#SBATCH --ntasks=1
# Request memory
#SBATCH --mem=20G
# Run for five minutes
#SBATCH --time=30:00:00
# Name job
#SBATCH --job-name=SbatchJob
# Name output file
#SBATCH --output=%x_%j.out

# change to the directory where you submitted this script
cd ${SLURM_SUBMIT_DIR}

# Executable section: echoing some Slurm data
echo "Starting sbatch script myscript.sh at:`date`"

cd /users/c/p/cpetak/WGS/geva

mychr=$1
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

read start stop < <(bcftools view -H ${input_vcf}_processed.recode.vcf.gz | awk 'NR==1 {first=$2} {last=$2} END {print first, last}')

python make_guide.py $start $mychr $stop

echo "Done with everything"
