#!/bin/bash

#Qualimap Counts_NJ
#qualimap: calculate the number of reads that mapp to a specific location
#run from pflaphy-ahcore/nadia/FRD3/

#$ -t 1-24
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N qualicounts_NJ
#$ -l vf=8G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 2


file=$(cat names.txt | head -n $SGE_TASK_ID | tail -n 1)

##qualimap <tool_name> <tool_options>
##qualimap comp-counts [-a <arg>] -bam <arg> -gtf <arg> [-id <arg>] [-out <arg>] [-p <arg>] [-pe] [-s <arg>] [-type <arg>]
/prj/pflaphy-ahcore/nadia/SOFTWARE/qualimap_v2.2.1/qualimap comp-counts -a proportional -bam /prj/pflaphy-ahcore/nadia/FRD3/Error_Cor/$file/$file'_outfile.bam' -gtf /prj/pflaphy-ahcore/nadia/GENOME/AHAL_GEM/annotation_Ahal_gem_onlyexons_t1_id.gtf -out /prj/pflaphy-ahcore/nadia/FRD3/Counts/Qualimap$file.txt -pe -p strand-specific-reverse
