#!/usr/bin/env bash

#BSUB -J Fastqc[1-2]
#BSUB -o logs/fastqc_%J.out
#BSUB -e logs/fastqc_%J.err
#BSUB -R "select[mem>16] rusage[mem=16] " 
#BSUB -q rna

module load fastqc

SAMPLES=(
KD_1
KD_2
)

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

r1_file=$(find ./ -regextype posix-extended -regex ".*/${sample}_S.*_L00.*_R1_001.fastq.gz")
r2_file=$(find ./ -regextype posix-extended -regex ".*/${sample}_S.*_L00.*_R2_001.fastq.gz")

echo $sample
echo $r1_file 
echo $r2_file

fastqc $r1_file $r2_file --outdir fastqc_test