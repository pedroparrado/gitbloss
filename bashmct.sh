#!/bin/bash --login

#SBATCH -N 1-1   
#SBATCH --tasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -t 44:00:00

module load /app/modules/languages/python/2.7.11

for ((i=3; i<=11; i+=4))
do
    python surfmct.py $i 1000 &
done
 

wait
