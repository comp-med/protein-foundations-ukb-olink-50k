#!/bin/sh

## extract proxies

#SBATCH --partition=compute
#SBATCH --job-name=feature_selection_variance
#SBATCH --account=sc-users
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-8757%50
#SBATCH --output=%x-%A-%2a.out

## change directory
cd <path>

## Use Array Index to select features
echo "Job ID: $SLURM_ARRAY_TASK_ID"
olink="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' input/variance.decomp.pop.proteins)"
pop="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $2}' input/variance.decomp.pop.proteins)"

echo "Node ID: $SLURM_NODELIST"

## Do some logging
echo "Protein: ${olink} and Population: ${pop}"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Starting script at: $date"

## This is the container to be used
R_CONTAINER='<path>/all_inclusive_rstudio_4.3.2.sif'

# This is the script that is executed
# Get with rstudioapi::getSourceEditorContext()$path
R_SCRIPT='<path>/03_feature_selection_updated.R'

# Enter all directories you need, simply in a comma-separated list
BIND_DIR="<path>"

## The container 
singularity exec \
--bind $BIND_DIR \
$R_CONTAINER Rscript $R_SCRIPT ${olink} ${pop}

## Some more logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Finishing script at: $date"
echo "Done!"
