#!/bin/bash
#SBATCH -t 0:30:00
#SBATCH -A berglandlab
#SBATCH -c 1
#SBATCH -p standard
#SBATCH --mem=100G
#SBATCH -N 1
#SBATCH -o lookup_mpl_%a.out

module load gcc/11.4.0
module load openmpi/4.1.4
module load R/4.3.1
p=$(grep '750000' *arr_${SLURM_ARRAY_TASK_ID}.out | cut -d = -f 2 | cut -d , -f 1)
Rscript /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/optim_precompute_mpl.R norm_sfs_${SLURM_ARRAY_TASK_ID}_20.txt 20 $p 0.2