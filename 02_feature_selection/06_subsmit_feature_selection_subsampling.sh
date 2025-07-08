#!/bin/sh

## extract proxies

#SBATCH --partition=compute
#SBATCH --job-name=feature_selection_variance_subsampling
#SBATCH --account=sc-users
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-3355%50
#SBATCH --output=%x-%A-%2a.out

## change directory
cd <path>
  
## Use Array Index to select features
echo "Job ID: $SLURM_ARRAY_TASK_ID"
olink="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var+19999 {print $1}' input/EUR.subsample.ancestry.txt)"
pop="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var+19999 {print $2}' input/EUR.subsample.ancestry.txt)"
grp="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var+19999 {print $3}' input/EUR.subsample.ancestry.txt)"

echo "Node ID: $SLURM_NODELIST"

## Do some logging
echo "Protein: ${olink} and Population: ${pop} using subsample ${grp}"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Starting script at: $date"

## This is the container to be used
R_CONTAINER='<path>/all_inclusive_rstudio_4.3.2.sif'

# This is the script that is executed
# Get with rstudioapi::getSourceEditorContext()$path
R_SCRIPT='<path>/02_feature_selection/scripts/07_feature_selection_updated_subsampling.R'

# Enter all directories you need, simply in a comma-separated list
BIND_DIR="<path>"

## The container 
singularity exec \
--bind $BIND_DIR \
$R_CONTAINER Rscript $R_SCRIPT ${olink} ${pop} ${grp}

## Some more logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Finishing script at: $date"
echo "Done!"
