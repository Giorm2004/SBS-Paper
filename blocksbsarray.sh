#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH -A berglandlab_standard
#SBATCH -c 1
#SBATCH -p standard
#SBATCH -N 1
module load gcc/11.4.0
module load openmpi/4.1.4
module load python/3.11.4

a=1
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ../test_control.txt)
for i in $OPTS; do declare "opt$a=$i"; ((a++)); done
echo $opt1
echo $opt2
echo $opt3
echo $opt4

slim=$'/home/nzx3cc/sbs_paper/apps/build/slim'
$slim -t -m -d ss=1.0 -d ds=0.6 -d sw=1 -d dw=0.6 -d r=1e-6 -d WinterFreq=0.5 -d len=20 -d N=1e4 -d name="${SLURM_ARRAY_TASK_ID}"  ../sbs.slim
python ../recap.py ./"${SLURM_ARRAY_TASK_ID}"_sbs_decap.trees 1e4 1e-6 1e-6 "${SLURM_ARRAY_TASK_ID}"_sbs
