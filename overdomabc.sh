#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH -c 40
#SBATCH -N 1
#SBATCH -A berglandlab
#SBATCH -p standard
#SBATCH --mem=100G
module load miniforge
conda activate testenv
#python abctest.py "blockoverdom/WFsim_overdomrecap_84.trees" 20000 30
python /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/abctest.py "./${SLURM_ARRAY_TASK_ID}"_overdom_recap.trees 10000 2 "${SLURM_ARRAY_TASK_ID}" 1e-6
