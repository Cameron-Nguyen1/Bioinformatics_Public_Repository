#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=80g
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH -o stdout.txt
#SBATCH -e error.txt
export R_LIBS_USER='/work/users/c/a/came/R_LIBS'
module load r/4.2.2
Rscript Sc_azimuth.R