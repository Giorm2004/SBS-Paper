#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH -A berglandlab
#SBATCH -c 1
#SBATCH -p standard
#SBATCH -N 1
#SBATCH --mem=9G
#SBATCH -o overdomarr_%a.out
module load gcc/11.4.0
module load openmpi/4.1.4
module load python/3.11.4

a=1
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/big_over_control.txt)
for i in $OPTS; do declare "opt$a=$i"; ((a++)); done
echo $opt1
echo $opt2
echo $opt3
name="${opt2}_${opt3}"
cd od_sims_$name

slim=$'/home/nzx3cc/sbs_paper/apps/build/slim'
$slim -t -m -d ss=-0.1 -d d=$opt3 -d sw=0 -d K=$opt2 -d r=1e-8 -d WinterFreq=0.5 -d len=20 -d N=1e4 -d name="${opt1}"  /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/overdom.slim
python /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/recap.py ./"${opt1}"_decap.trees 1e4 1e-8 1e-8 "${opt1}"_