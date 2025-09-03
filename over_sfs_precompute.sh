#!/bin/bash
#SBATCH -t 0:15:00
#SBATCH -A berglandlab
#SBATCH -c 1
#SBATCH -p standard
#SBATCH -N 1
#SBATCH --mem=100G
REAL_ID=$(($1 + $SLURM_ARRAY_TASK_ID))
echo $REAL_ID
a=1
OPTS=$(sed -n "${REAL_ID}"p /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/over_sfs_control.txt)
for i in $OPTS; do declare "opt$a=$i"; ((a++)); done
echo $opt1
echo $opt2
echo $opt3
echo $opt4

module load gcc/11.4.0
module load openmpi/4.1.4
module load R/4.3.1
Rscript /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/overdom_precompute.R $opt1 $opt2 $opt3 $opt4