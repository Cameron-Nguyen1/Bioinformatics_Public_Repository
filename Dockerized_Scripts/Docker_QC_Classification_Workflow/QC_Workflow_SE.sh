#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=80g
#SBATCH -n 12
#SBATCH -t 48:00:00

#Load libraries for processing data
module load kraken
module load gcc/11.2.0
module load bowtie2
module load samtools
module load fastp

#Step 1. fastp QC
#for r1 in *R1*.fastq #PE_DATA
for r1 in *.fastq #SE_DATA
do
        #r2=$(echo $r1 | sed -e 's/R1/R2/g')
        #fastp --dedup -c -w 12 -i $r1 -I $r2 -o ${r1::-6}".trim.fastq" -O ${r2::-6}".trim.fastq" #PE_DATA
        fastp --dedup -c -w 12 -i $r1 -o ${r1::-6}".trim.fastq" #SE_DATA
done

#Step 2. Bowtie2 Alignment
#for R1 in *R1*.trim.fastq #PE
for r1 in *.trim.fastq
do
    #R2=$(echo $R1 | sed -e 's/R1/R2/g') #PE
    out=${r1::-11}
    bowtie2 --phred33 -q -p 12 -x /work/users/c/a/came/Virs/Reads/bt2_genome/SARS2 -U $r1 -S $out".output.sam" --un $out".unmapped.fastq"
    samtools view -S -b $out".output.sam" > $out".output.bam"
    samtools sort $out".output.bam" -o $out".sorted.bam"
    samtools index $out".sorted.bam"
    continue
done

for var in *.unmapped.fastq
do
        kraken2 --use-names --threads 12 --db /nas/longleaf/apps/kraken2/2.1.2/KRAKEN2-HS-FUNGI-MM-DB --report ${var::-40}".report" $var > ${var::-40}".kraken"
done