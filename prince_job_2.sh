#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=10:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=c-2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ca2356@nyu.edu
#SBATCH --output=slurm_%j.out

methods=($(cat methods.txt))
method=${methods[1]}
#./parallel_execute.sh L10_$method
./parallel_execute.sh L10plus_$method
