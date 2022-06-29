#!/bin/bash

#Error Correction_NJ
#script_comex
#run from pflaphy-ahore/nadia/VPA/ErrCorr

#$ -t 1-24
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N errorCorrect_NJ
#$ -l vf=40G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat names.txt | head -n $SGE_TASK_ID | tail -n 1)

direct=$file
mkdir $direct
echo "Directory for $direct created"

InputSam=/prj/pflaphy-ahcore/nadia/FRD3/Mapping/$direct.sam
if [ -f $InputSam ];
	then
		echo "File $InputSam exists."
	else
		echo "File $InputSam does not exist."
		exit
fi

IFS=\.

arr=($InputSam)


input=$direct
echo $input

IFS=$OIFS

referenceGenome=/prj/pflaphy-ahcore/nadia/GENOME/AHAL_GEM/Halleri.renamed.fasta

if [ -f $referenceGenome ];
	then
		echo "File $referenceGenome exists."
	else
		echo "File $referenceGenome does not exist."
		exit
fi

python ./toprintend1.py /prj/pflaphy-ahcore/nadia/FRD3/Mapping/$input'.sam' ./$direct/$input'_end.sam'

echo "End created"

python ./Selectnonrepeated1.py ./$direct/$input'_end.sam' ./$direct/$input'_hits_nonrepeated.sam'

echo "non-repeated hits finished"

python ./removeEnd1.py ./$direct/$input'_hits_nonrepeated.sam' >./$direct/$input'_output_final.sam'

echo "End removed"

samtools view -bT $referenceGenome ./$direct/$input'_output_final.sam' >./$direct/$input'_outfile.bam'

echo "All tophat-errors corrected for $input."

rm ./$direct/$input'_end.sam'
rm ./$direct/$input'_hits_nonrepeated.sam'

echo "Temporary files have been deleted"


