#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH -c 40
#SBATCH -N 1
#SBATCH -A berglandlab_standard
#SBATCH -p standard
#SBATCH --mem=50G
module load miniforge
conda activate testenv
#python abctest.py "blockoverdom/WFsim_overdomrecap_84.trees" 20000 30
python ../abctest.py "./${SLURM_ARRAY_TASK_ID}"_overdom_recap.trees 20000 30 "${SLURM_ARRAY_TASK_ID}" 40
