#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --time=10:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=coupled-1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ca2356@nyu.edu
#SBATCH --output=slurm_%j.out

#BAOAB
#memory_BAOAB
#middle_IsoK_NHC
#middle_NHC
#middle_NHL
./parallel_execute.sh middle_SINR
./parallel_execute.sh middle_SubK_NHC
./parallel_execute.sh xi_respa_NHC
./parallel_execute.sh xo_respa_NHC
