#!/bin/sh

## submit imputation

#SBATCH --partition=compute
#SBATCH --job-name=impute_data
#SBATCH --account=sc-users
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=110
#SBATCH --mem-per-cpu=3G
#SBATCH --array=1
#SBATCH --output=%x-%A-%2a.out

## change directory
cd <path>

## Use Array Index to select features
echo "Job ID: $SLURM_ARRAY_TASK_ID"

## Do some logging
echo "Perform impuation"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Starting script at: $date"

## This is the container to be used
R_CONTAINER='<path>/all_inclusive_rstudio_4.3.2.sif'

# This is the script that is executed
# Get with rstudioapi::getSourceEditorContext()$path
R_SCRIPT='<path>/03_perform_imputation.R'

# Enter all directories you need, simply in a comma-separated list
BIND_DIR="<path>"

## The container 
singularity exec \
--bind $BIND_DIR \
$R_CONTAINER Rscript $R_SCRIPT

## Some more logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Finishing script at: $date"
echo "Done!"
