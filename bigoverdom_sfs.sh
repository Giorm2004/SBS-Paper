#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH -A berglandlab
#SBATCH -c 1
#SBATCH -p standard
#SBATCH -N 1
#SBATCH --mem=9G
#SBATCH -o sfs_%a.out
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

python /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/get_sfs.py "${opt1}"__recap.trees 20 100000 200 ${opt1} 20000 0 1000
python /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/get_sfs.py "${opt1}"__recap.trees 20 100000 200 ${opt1} 20000 20000 1000
python /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/get_sfs.py "${opt1}"__recap.trees 20 100000 200 ${opt1} 20000 40000 1000