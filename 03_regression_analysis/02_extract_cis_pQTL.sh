#!/bin/sh

## script to extract top cis-pQTL across different ancestries from UKB-PPP stats

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
#SBATCH --mem-per-cpu=5G

#! how many cpus per task
#SBATCH --cpus-per-task=10

#! no other jobs can be run on the same node
##SBATCH --exclusive

#! run as job array
#SBATCH --array=53

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=slurm-%x-%j.out

## change directory
cd <path>
  
## get the exposure
echo "Job ID: $SLURM_ARRAY_TASK_ID"
olink="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var+1 {print $3}' input/cis.pQTL.query.ancestries.txt)"
chr="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var+1 {print $7}' input/cis.pQTL.query.ancestries.txt)"
start="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var+1 {print $8}' input/cis.pQTL.query.ancestries.txt)"
end="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var+1 {print $9}' input/cis.pQTL.query.ancestries.txt)"

echo "Node ID: $SLURM_NODELIST"

echo "Protein: ${olink} | chromosome : ${chr} | start : ${start} | end : ${end}"

## --> extract statistics <-- ##

for region in EUR AFR CSA; do

  echo "Processing region: $region"
  
  ## get folder
  folder="$(ls <path>/UKB_PPP/${region}/ | grep -v 'tar' | grep "^${olink}_")"
  ## print
  echo "$folder"
  
  ## get file
  echo "chr${chr}_"
  file="$(ls <path>/UKB_PPP/${region}/$folder | grep "chr${chr}_")"
  ## print
  echo "$file"
  
  ## extract the needed information; drop low-freq and poor quality variants
  zcat "<path>/UKB_PPP/${region}/${folder}/${file}" | \
  awk -v low=${start} -v upp=${end} '{if(($2+0 >= low-5e5 && $2+0 <= upp+5e5 && $6+0 >= 0.005 && $6+0 <= 0.995 && $7+0 >= 0.4) || NR == 1) print $0}' > input/${region}.${olink}.cis.region.txt

done
