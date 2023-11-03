#!/bin/bash
#SBATCH -p general
#SBATCH -N X
#SBATCH --mem=Y
#SBATCH -n 12
#SBATCH -t 48:00:00

#Load libraries for processing data
module load kraken
module load gcc/11.2.0
module load bowtie2
module load samtools
module load fastp

#Step 1. fastp QC

for r1 in *.fastq #SE_DATA
do
        fastp --dedup -c -w 12 -i $r1 -o ${r1::-6}".trim.fastq" #SE_DATA
done

#Step 2. Bowtie2 Alignment

for r1 in *.trim.fastq
do
    out=${r1::-11}
    bowtie2 --phred33 -q -p 12 -x /path/to/bt2_genome/SARS2 -U $r1 -S $out".output.sam" --un $out".unmapped.fastq"
    samtools view -S -b $out".output.sam" > $out".output.bam"
    samtools sort $out".output.bam" -o $out".sorted.bam"
    samtools index $out".sorted.bam"
    continue
done

#Step 3. Kraken2 classification 
for var in *.unmapped.fastq
do
        kraken2 --use-names --threads 12 --db /path/to/kraken2_db/db --report ${var::-40}".report" $var > ${var::-40}".kraken"
done