#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=16g
#SBATCH -n 2
#SBATCH -t 1-
module load r/4.2.2
#Rscript Kali_Dependencies.R #Comment this out if you've already ran the Dependency installation script
base=$(pwd)
export R_LIBS=/example/path/R_LIBS:$R_LIBS 

Rscript ks6.R --studydesign Study_Design.txt --database EnsDb.Mmusculus.v79 \
--expression "infection=Group1-Group2" \
--ex_uh "(Group2_1+Group2_2+Group2_3+Group2_4+Group2_5)/5" \
--ex_eh "(Group1_1+Group1_2+Group1_3+Group1_4+Group1_5)/5" \
--species "Mus musculus" \
--out ${base}"/Output"