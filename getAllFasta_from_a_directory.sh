#!/bin/bash
printf "%s\n" "Let's get fasta!" ""
SPIN='-\|/'
ls | grep ".bed" > filename.tmp
numfile=$(cat filename.tmp | wc -l | sed 's/ //g')

if [[ $numfile -eq 0 ]]; then
	printf "Error: The current directory does not present bed file!\n\n"
elif [[ $numfile -eq 1 ]]; then
	printf "The current directory presents 1 bed file\n\n"
else
	printf "The current directory presents ${numfile} bed files\n\n"
fi
for (( k = 1; k < ${numfile} +1; k++ )); do
	filename=$(sed -n $k'p' filename.tmp)
	bedtools getfasta -fi /Users/boulanger/Desktop/Resultats/ChIP/Chip_seq_analysis/Annotation_genome/GRCh37_Hg19/GRCh37_20190212_mainChr_chrFormat.fna -bed ${filename} -fo ${filename%%.*}.fa
done & PID=$!
i=0 &
while kill -0 $PID 2>/dev/null; do
	i=$(( (i+1) %4 ))
	printf "\rFasta recovery in process ${SPIN:$i:1} "
	sleep .1
done
printf "\nAll sequences contained in bed files had been recovered in fasta file!\n"
rm filename.tmp
exit
