#!/bin/bash
#SBATCH -N 1                # number of nodes
#SBATCH -n 16               # number of cores
#SBATCH -t 48:00:00
#SBATCH --job-name=inlets
#SBATCH --mem 64g
#SBATCH -o stdout.txt
#SBATCH -e error.txt

module load cellranger

cellranger multi --id=REP_C --localmem=64 --localcores=16 \
--csv=config_inlet.csv