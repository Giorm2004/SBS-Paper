#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH -A berglandlab
#SBATCH -c 1
#SBATCH -p standard
#SBATCH -N 1
#SBATCH --mem=9G
#SBATCH -o mpl_%a.out
module load gcc/11.4.0
module load openmpi/4.1.4
module load R/4.3.1


a=1
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/big_sbs_control.txt)
for i in $OPTS; do declare "opt$a=$i"; ((a++)); done
echo $opt1
echo $opt2
echo $opt3
name="${opt2}_${opt3}"


END=$(($opt2 * 10000))
T1=$(($END - 5))
T2=$(($END -10))

p1=$(grep ${END}: *arr_${SLURM_ARRAY_TASK_ID}.out | cut -d \; -f 2| cut -d = -f 2 | cut -d , -f 1)
p2=$(grep ${T1}: *arr_${SLURM_ARRAY_TASK_ID}.out | cut -d = -f 2 | cut -d , -f 1)
p3=$(grep ${T2}: *arr_${SLURM_ARRAY_TASK_ID}.out | cut -d = -f 2 | cut -d , -f 1)

echo sbs_sims_$name
cd sbs_sims_$name

echo $p1
echo $p2
echo $p3
Rscript /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/optim_mpl.R norm_sfs_${opt1}_20_40000.txt 20 $p1 0.15
Rscript /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/optim_mpl.R norm_sfs_${opt1}_20_20000.txt 20 $p2 0.15
Rscript /scratch/nzx3cc/nzx3cc/sbs_paper/scripts/optim_mpl.R norm_sfs_${opt1}_20_0.txt 20 $p2 0.15