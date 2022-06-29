#!/bin/bash

#Hisat2
#alignment_rnaseq_NJ
#run from pflaphy-ahcore/nadia/FRD3

#$ -t 1-24
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N map_rna_NJ
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat names.txt | head -n $SGE_TASK_ID | tail -n 1)

## perl /path/to/hisat2 -x /path/to/IndexName(without extension) -1 /path/to/lib1.for.paired.fq.gz -2 /path/to/lib1.rev.paired.fq.gz -S /path/to/lib1.aligned.sam --summary-file /path/to/lib1.summary.txt
## hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
/prj/pflaphy-ahcore/nadia/SOFTWARE/hisat2-2.1.0/hisat2 --no-unal --rna-strandness RF -q -x /prj/pflaphy-ahcore/nadia/GENOME/AHAL_GEM/Index_AHalleri/Halleri.renamed.index -1 /prj/pflaphy-ahcore/nadia/FRD3/Trimmed/$file.1.cutadapt.fastq.gz -2 /prj/pflaphy-ahcore/nadia/FRD3/Trimmed/$file.2.cutadapt.fastq.gz -S /prj/pflaphy-ahcore/nadia/FRD3/Mapping/$file.sam --summary-file /prj/pflaphy-ahcore/nadia/FRD3/Mapping/Summary/$file-summary.txt --max-intronlen 20000 -L 10 --pen-noncansplice 16 --mp 2,0 --score-min L,0,-0.3 --rdg 3,1 --rfg 3,1 --pen-canintronlen G,-0.5,0.1 2> /prj/pflaphy-ahcore/nadia/FRD3/Mapping/$file.align.err

echo "Alignment done"

