#!/bin/sh

## obtain SNPs

#SBATCH --partition=compute
#SBATCH --job-name=obtain_snps
#SBATCH --account=sc-users
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-132%5
#SBATCH --output=%x-%A-%2a.out

## change directory
cd <path>
  
## Do some logging
echo "Perform impuation"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Starting script at: $date"

## export location of files
export dir=<path>
export out=<path>

## Use Array Index to select features
echo "Job ID: $SLURM_ARRAY_TASK_ID"
chr="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $2}' files.obtain.snps.txt)"
chunk="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $3}' files.obtain.snps.txt)"

## what is done
echo "Obtain SNPs from chromosome ${chr} using chunk ${chunk}"

## create subset bgen file
<path>/bgenix \
-g ${dir}/ukb22828_c${chr}_b0_v3.bgen \
-incl-rsids snps.${chr}.${chunk}.txt > ${out}/tmp.${chr}.${chunk}.bgen

## create dosage file
<path>/qctool_v2.0.7 \
-g ${out}/tmp.${chr}.${chunk}.bgen \
-s ${dir}/ukb22828_c${chr}_b0_v3.sample  \
-og - \
-ofiletype dosage > ${out}/tmp.${chr}.${chunk}.dosage

cd ${out}

## get SNP info
cut -f 1-6 -d ' ' tmp.${chr}.${chunk}.dosage > snp.${chr}.${chunk}.info

## transpose dosage matrix (https://stackoverflow.com/questions/1729824/an-efficient-way-to-transpose-a-file-in-bash)
cut -f 1-6 -d ' ' --complement tmp.${chr}.${chunk}.dosage | awk '
{
  for (i=1; i<=NF; i++)  {
      a[NR,i] = $i
  }
}
NF>p { p = NF }
END {
  for(j=1; j<=p; j++) {
      str=a[1,j]
      for(i=2; i<=NR; i++){
          str=str" "a[i,j];
      }
      print str
  }
}' - > snp.${chr}.${chunk}.dosage.transpose


## delete files not longer needed
rm tmp.${chr}.${chunk}.dosage
rm tmp.${chr}.${chunk}.bgen

## Some more logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Finishing script at: $date"
echo "Done!"
