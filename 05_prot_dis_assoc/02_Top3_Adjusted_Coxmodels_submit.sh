#!/bin/sh

## run simple Cox models with proteins for prot variance project

#SBATCH --partition=compute
#SBATCH --job-name=Top3_Adj_Cox_Protvar
#SBATCH --account=sc-users
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-443%50
#SBATCH --output=/sc-projects/sc-proj-computational-medicine/people/Kamil/projects/03_UKB_Protein_Variance/01_Cox_prediction/slurm_logs/%x-%A-%2a.out

## change directory
cd /sc-projects/sc-proj-computational-medicine/people/Kamil/projects/03_UKB_Protein_Variance/01_Cox_prediction

## Use Array Index to select features
echo "Job ID: $SLURM_ARRAY_TASK_ID"
phe="$(awk -v var=$SLURM_ARRAY_TASK_ID -F '\t' 'NR == var {print $1}' input/input_phecodes.txt)"
cohort="$(awk -v var=$SLURM_ARRAY_TASK_ID -F '\t' 'NR == var {print $2}' input/input_phecodes.txt)"

## Do some logging
echo "Running prediction for ${phe} in ${cohort}"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Starting script at: $date"

## This is the container to be used
R_CONTAINER='/sc-projects/sc-proj-computational-medicine/programs/all-inclusive-rstudio-apptainer/sif/all_inclusive_rstudio_4.3.2.sif'

# This is the script that is executed
# Get with rstudioapi::getSourceEditorContext()$path
R_SCRIPT='/sc-projects/sc-proj-computational-medicine/people/Kamil/projects/03_UKB_Protein_Variance/01_Cox_prediction/scripts/02_Top3_Adjusted_Coxmodels.R'

# Enter all directories you need, simply in a comma-separated list
BIND_DIR="/sc-projects/sc-proj-computational-medicine/,/sc-resources/"

## The container 
singularity exec \
--bind $BIND_DIR \
$R_CONTAINER Rscript $R_SCRIPT ${phe} ${cohort}

## Some more logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Finishing script at: $date"
echo "Done!"