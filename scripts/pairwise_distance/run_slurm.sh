#!/bin/bash

#Submit this script with:sbatch thefilename

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=1G   # memory per CPU core
#SBATCH -J "snakemaster"   # job name


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load julia/1.8.2

snakemake --use-conda --cores=80 --cluster 'sbatch -t 3000 --ntasks={threads} --cores={threads} --mem-per-cpu={resources.mem_mb} --nodes=1' -j 80 --scheduler greedy
