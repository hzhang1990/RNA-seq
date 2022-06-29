#!/bin/bash

#Quality check_before trimming
#FastQC
#run from /prj/pflaphy-ahcore/nadia/FRD3 

#$ -t 1-48
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N qc_bt_NJ
#$ -l vf=2G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12


file=$(cat /prj/pflaphy-ahcore/nadia/FRD3/list.txt | head -n $SGE_TASK_ID | tail -n 1)
/prj/pflaphy-ahcore/nadia/SOFTWARE/FastQC/fastqc -t 4 /prj/pflaphy-ahcore/nadia/FRD3/Raw_data/$file -o /prj/pflaphy-ahcore/nadia/FRD3/FastQC_BT/
#md5sum /prj/pflaphy-ahcore/nadia/FRD3/Raw_data/$file > /prj/pflaphy-ahcore/nadia/FRD3/FastQC_BT/md5sum/$file.md5

