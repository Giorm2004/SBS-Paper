#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH -A berglandlab
#SBATCH -c 1
#SBATCH -p standard
#SBATCH -N 1
#SBATCH -o neutarr_%a.out
module load gcc/11.4.0
module load openmpi/4.1.4
module load python/3.11.4

#slim=$'/home/nzx3cc/sbs_paper/apps/build/slim'
#$slim -t -m -d r=1e-8 -d N=1e4 -d name="${SLURM_ARRAY_TASK_ID}"  /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/neut.slim
python /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/recap.py ./"${SLURM_ARRAY_TASK_ID}"neut_decap.trees 1e4 1e-8 1e-8 "${SLURM_ARRAY_TASK_ID}"_neut