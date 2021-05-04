#!/usr/bin/env bash

#BSUB -J trim_bowtie[1-2]
#BSUB -o logs/trim_bowtie_%J.out
#BSUB -e logs/trim_bowtie_%J.err
#BSUB -R "select[mem>16] rusage[mem=16] " 
#BSUB -q rna
#BSUB -n 6

module load trimgalore/0.4.5
module load bowtie2/2.3.2

SAMPLES=(
Input_1_S1
Input_2_S2
)

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
indexes=indexes/mm10
trim_dir=trim
aligned_dir=aligned

r1_file=${sample}_L004_R1_001.fastq.gz
r2_file=${sample}_L004_R2_001.fastq.gz
r1_trim=$trim_dir/${sample}_L004_R1_001_val_1.fq.gz
r2_trim=$trim_dir/${sample}_L004_R2_001_val_2.fq.gz
sam_output=$aligned_dir/${sample}.sam

echo $sample

mkdir $aligned_dir
mkdir $trim_dir
trim_galore --paired -o $trim_dir $r1_file $r2_file
bowtie2 --very-sensitive -p 6 -x $indexes -1 $r1_trim -2 $r2_trim -S $sam_output