#!/bin/sh

## extract proxy variants

#! select partition
#SBATCH --partition=compute

#! select account
#SBATCH --account=sc-users

#! Specify required run time
#SBATCH --time=5:00:00

#! how many nodes
#SBATCH --nodes=1

#! how many tasks
#SBATCH --ntasks=1

#! how much memory per cpu
#SBATCH --mem-per-cpu=20G

#! how many cpus per task
#SBATCH --cpus-per-task=20

#! no other jobs can be run on the same node
##SBATCH --exclusive

#! run as job array
#SBATCH --array=1-23%3

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=slurm-%x-%j.out

## change directory
cd <path>

## get the exposure
echo "Job ID: $SLURM_ARRAY_TASK_ID"

## export location of genotype files to be used
export dir=<path>
## directory of variant inclusion files
export var=<path>

## run based on which chromosome is choosen
if [[ $SLURM_ARRAY_TASK_ID -eq 23 ]]; then

## define X-chromsome
export chr="X"

## run plink to obtain matrix of SNPs in LD
<path>/plink2_20241025 \
--pfile ${dir}/plink2_chr${chr} \
--r2-unphased \
--keep <path>/random_sample/UKBB.unrelated.random.30k.2024-06-13.sample \
--ld-window-kb 3000 \
--ld-window-r2 0.5 \
--threads 20 \
--ld-snp-list snp.list.23.txt \
--out ld.proxies.23

else
  
  ## any chromosome
  export chr=$SLURM_ARRAY_TASK_ID

## run plink to obtain matrix of SNPs in LD
<path>/plink2_20241025 \
--pfile ${dir}/plink2_chr${chr} \
--r2-unphased \
--keep <path>/random_sample/UKBB.unrelated.random.30k.2024-06-13.sample \
--ld-window-kb 3000 \
--ld-window-r2 0.5 \
--threads 20 \
--ld-snp-list snp.list.${chr}.txt \
--out ld.proxies.${chr}

fi
