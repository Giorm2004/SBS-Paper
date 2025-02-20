#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH -A berglandlab
#SBATCH -c 1
#SBATCH -p standard
#SBATCH -N 1
#SBATCH --mem=9G
module load gcc/11.4.0
module load openmpi/4.1.4
module load python/3.11.4

a=1
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/simcontrol.txt)
for i in $OPTS; do declare "opt$a=$i"; ((a++)); done
echo $opt1
echo $opt2
echo $opt3
echo $opt4

slim=$'/home/nzx3cc/sbs_paper/apps/build/slim'
$slim -t -m -d ss=$opt1 -d d=1.5 -d sw=0.5 -d r=1e-7 -d WinterFreq=0.5 -d len=20 -d N=1e4 -d name="${SLURM_ARRAY_TASK_ID}"  /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/overdom.slim
python /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/recap.py ./"${SLURM_ARRAY_TASK_ID}"overdom_decap.trees 1e4 1e-6 1e-7 "${SLURM_ARRAY_TASK_ID}"_overdom
