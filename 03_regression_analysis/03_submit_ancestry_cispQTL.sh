#!/bin/sh

## run cross-ancestry evaluation of cispQTLs

#SBATCH --partition=compute
#SBATCH --job-name=cispQTL_ancestry
#SBATCH --account=sc-users
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=15G
#SBATCH --array=3-137%10
#SBATCH --output=%x-%A-%2a.out

## change directory
cd <path>/03_regression_analysis

## Use Array Index to select features
echo "Job ID: $SLURM_ARRAY_TASK_ID"
olink="$(awk -v var=$SLURM_ARRAY_TASK_ID -F '\t' 'NR == var+1 {print $3}' input/cis.pQTL.query.ancestries.txt)"

## Do some logging
echo "Running extraction for ${olink}"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Starting script at: $date"

## This is the container to be used
R_CONTAINER='<path>/all_inclusive_rstudio_4.3.2.sif'

# This is the script that is executed
# Get with rstudioapi::getSourceEditorContext()$path
R_SCRIPT='<path>/04_cis_pQTL_testing.R'

# Enter all directories you need, simply in a comma-separated list
BIND_DIR="<path>"

## The container 
singularity exec \
--bind $BIND_DIR \
$R_CONTAINER Rscript $R_SCRIPT ${olink}

## Some more logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Finishing script at: $date"
echo "Done!"
