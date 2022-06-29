#!/bin/bash

#Trimming
#cutadapt 2.1
#run from pflaphy-ahcoe/nadia/FRD3/

#$ -t 1-24
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N cutadapt_frd3rna_NJ
#$ -l vf=2G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi


file=$(cat names.txt | head -n $SGE_TASK_ID | tail -n 1)

## For paired-end reads: cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq

/vol/python/bin/cutadapt -e 0.15 -O 4 -m 120 --nextseq-trim=20 -q 20 -a 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC' -A 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT' -a "A{100}" -n 2 -o ./Trimmed/$file.1.cutadapt.fastq.gz -p ./Trimmed/$file.2.cutadapt.fastq.gz ./Raw_data/$file/$file'_1.fq.gz' ./Raw_data/$file/$file'_2.fq.gz' 1> ./Trimmed/reports/$file.cutadapt 2> ./Trimmed/reports/$file.cutadapt.err

#chmod 755 cutadapt.sh
#nohup ./cutadapt.sh &
#option -j number of cores for multithreading available together with pigz installed!
