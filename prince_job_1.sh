#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=10:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=c-1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ca2356@nyu.edu
#SBATCH --output=slurm_%j.out

./parallel_execute.sh BAOAB
./parallel_execute.sh memory_BAOAB
./parallel_execute.sh middle_IsoK_NHC
#./parallel_execute.sh middle_NHC
#./parallel_execute.sh middle_NHL
#middle_SINR
#middle_SubK_NHC
#xi_respa_NHC
#xo_respa_NHC
