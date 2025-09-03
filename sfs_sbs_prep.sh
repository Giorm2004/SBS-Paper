#!/bin/bash
#SBATCH -t 0:20:00
#SBATCH -A berglandlab
#SBATCH -c 1
#SBATCH -p standard
#SBATCH -N 1
#SBATCH --mem=100G
module load gcc/11.4.0
module load openmpi/4.1.4
module load python/3.11.4
python /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/get_sfs.py "${SLURM_ARRAY_TASK_ID}"_sbs_recap.trees $1 100000 200 ${SLURM_ARRAY_TASK_ID} 1000