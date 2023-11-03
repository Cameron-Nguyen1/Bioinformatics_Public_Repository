#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=80g
#SBATCH -n 12
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=came@unc.edu # send-to address
#SBATCH -o QC_stdout_count.txt
#SBATCH -e QC_stderr_count.txt

#Load libraries for processing data
module load kraken
module load gcc/11.2.0
module load bowtie2
module load samtools
module load fastp

#Step 1. fastp QC
for r1 in *R1*.fastq
do
        r2=$(echo $r1 | sed -e 's/R1/R2/g')
        fastp --dedup --dup_calc_accuracy 6 --detect_adapter_for_pe -5 -c -p -w 8 -i $r1 -I $r2 -o ${r1::-6}".trim.fastq" -O ${r2::-6}".trim.fastq" -j ${r1::-6}".fastp.json" -h ${r1::-6}".fastp.html"
done

#Step 2. Bowtie2 Alignment + Consensus + VCF
for R1 in *R1*.trim.fastq
do
        R2=$(echo $R1 | sed -e 's/R1/R2/g')
        out=${R1::-15}

        if [[ "$R1" == *"OC43"* ]]
        then
                bowtie2 --phred33 -q -p 12 -x /work/users/c/a/came/Virs/Reads/bt2_genome/OC43 -1 $R1 -2 $R2 -S $out".output.sam" --un-conc $out".unmapped.fastq"
                samtools view -S -b $out".output.sam" > $out".output.bam"
                samtools sort $out".output.bam" -o $out".sorted.bam"
                samtools index $out".sorted.bam"
                bcftools mpileup -O b -o $out".bcf" -f /work/users/c/a/came/Virs/Reads/bt2_genome/OC43.fasta $out".sorted.bam"
                bcftools call --ploidy 1 -m -v -o $out".vcf" $out".bcf" 
                vcfutils.pl varFilter $out".vcf" > $out"_final_variants.vcf"
                bgzip $out"_final_variants.vcf"
                bcftools index $out"_final_variants.vcf.gz"
                bcftools consensus -f /work/users/c/a/came/Virs/Reads/bt2_genome/OC43.fasta $out"_final_variants.vcf.gz" > $out"_consensus.fasta"
                continue

        elif [[ "$R1" == *"MA10"* ]]
        then
                bowtie2 --phred33 -q -p 12 -x /work/users/c/a/came/Virs/Reads/bt2_genome/MA10 -1 $R1 -2 $R2 -S $out".output.sam" --un-conc $out".unmapped.fastq"
                samtools view -S -b $out".output.sam" > $out".output.bam"
                samtools sort $out".output.bam" -o $out".sorted.bam"
                samtools index $out".sorted.bam"
                bcftools mpileup -O b -o $out".bcf" -f /work/users/c/a/came/Virs/Reads/bt2_genome/MA10.fasta $out".sorted.bam"
                bcftools call --ploidy 1 -m -v -o $out".vcf" $out".bcf" 
                vcfutils.pl varFilter $out".vcf" > $out"_final_variants.vcf"
                bgzip $out"_final_variants.vcf"
                bcftools index $out"_final_variants.vcf.gz"
                bcftools consensus -f /work/users/c/a/came/Virs/Reads/bt2_genome/MA10.fasta $out"_final_variants.vcf.gz" > $out"_consensus.fasta"
                continue

        fi
done

for var in *.unmapped.1.fastq
do
        var2=$(echo $var | sed -e 's/1.fastq/2.fastq/g')
        kraken2 --use-names --paired --threads 12 --db /nas/longleaf/apps/kraken2/2.1.2/KRAKEN2-HS-FUNGI-MM-DB --report ${var::-20}".report" $var $var2 > ${var::-20}".kraken"
done