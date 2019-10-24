#!/bin/bash
#Made by Mathias Boulanger - 2019/02/25
#
#USAGE: gfforbed_to_pourcentage_coverage_genome.sh Bed_orGff.file human_genome_version (hg19 | hg38)
#version 1
#use on gff3 or bed file structure

ARGS=2							#The script need 1 argument
NAMEPROG=$(basename ${0})		#Name of the programme
DATA=$1							#File in argument
EXTENSION1="gff"				#Extension file necessary to run this script
EXTENSION2="bed"				#Extension file necessary to run this script
SPIN='-\|/'						#Waiting characters
RED='\033[1;31m'
GREEN='\033[1;32m'
ORANGE='\033[0;33m'
NOCOLOR='\033[0m'

##Check the ability to work
if [[ $# -ne $ARGS ]]; then
    printf "%s/n" "Wrong number of arguments!" "Usage: ${NAMEPROG} target_file.gff Human_genome_version (hg19 | hg38)" "Usage: ${NAMEPROG} target_file.bed Human_genome_version (hg19 | hg38)"
    exit 1
fi
if [[ ! -f $DATA ]];then
	printf "%s\n" "${RED}Error:${NOCOLOR} file does not exit!" "Please use an existing file." "Usage: ${NAMEPROG} target_file.gff Human_genome_version (hg19 | hg38)" "Usage: ${NAMEPROG} target_file.bed Human_genome_version (hg19 | hg38)"
	exit 1
fi
if [[ ${DATA##*.} != $EXTENSION1 && ${DATA##*.} != $EXTENSION2 ]]; then
	printf "\n"
	printf "${ORANGE}WARNING:${NOCOLOR} The file extension sould be ." $EXTENSION1 "or ." $EXTENSION2
	printf "Make sure that your file presents an gff or bed structure." "\n"
fi
if [[ $(wc -l $DATA) = "0 ${DATA}" ]]; then
	printf "${RED}Error:${NOCOLOR} your file is empty!\n"
	exit 1
fi
#source NCBI assembly
if [[ ${2} == "hg19" ]]; then
	BASEHUMANGENOME="3101788170"
elif [[ ${2} == "hg38" ]]; then
	BASEHUMANGENOME="3099706404"
else
	printf "%s\n" "${RED}Error:${NOCOLOR} Wrong version of the human genome!" "Usage: ${NAMEPROG} target_file.gff Human_genome_version (hg19 | hg38)" "Usage: ${NAMEPROG} target_file.bed Human_genome_version (hg19 | hg38)"
	exit 1
fi

##Start to work
#for gff file
if [[ ${DATA##*.} == $EXTENSION1 ]]; then
	if [[ $(awk 'BEGIN{FS="\t"}{print NF}' $DATA | uniq ) -ne 9 ]]; then
		printf "%s\n" "" "ERROR: The structure of your file (${DATA}) does not present the classical structure of gff file" "Please try with an orther file"
		exit 1
	else
		awk 'BEGIN{FS="\t"; OFS="\t"}{print ($5-$4)+1}' $DATA | awk 'BEGIN{FS="\t"; OFS="\t"}{SUM += $1}END{print (SUM/"'$BASEHUMANGENOME'")*100}' > .pourcent.tmp & PID=$!
		i=0 &
		while kill -0 $PID 2>/dev/null; do
			i=$(( (i+1) %4 ))
			printf "\rCalculation in process ${SPIN:$i:1}"
		sleep .1
		done
		POURCENTFILE="$(cat .pourcent.tmp)"
		rm .pourcent.tmp
		printf "\n\nThe coverage of the regions in your GFF file represents "$POURCENTFILE"%% of the human genome GRCh37p13 (hg19).\n\n"
		exit 0
	fi
elif [[ $(awk 'BEGIN{FS="\t"}{print NF}' $DATA && uniq) -eq 9 ]]; then
	awk 'BEGIN{FS="\t"; OFS="\t"}{print ($5-$4)+1}' $DATA | awk 'BEGIN{FS="\t"; OFS="\t"}{SUM += $1} END{print (SUM/"'$BASEHUMANGENOME'")*100}' > .pourcent.tmp & PID=$!
	i=0 &
	while kill -0 $PID 2>/dev/null; do
		i=$(( (i+1) %4 ))
		printf "\rCalculation in process ${SPIN:$i:1}"
	sleep .1
	done
	POURCENTFILE="$(cat .pourcent.tmp)"
	rm .pourcent.tmp
	printf "\n\nThe coverage of the regions in your GFF file represents "$POURCENTFILE"%% of the human genome GRCh37p13 (hg19).\n\n"
	exit 0
fi

#for bed file
if [[ ${DATA##*.} == $EXTENSION2 ]]; then
	if [[ $(awk 'BEGIN{FS="\t"}{print NF}' $DATA | uniq) -ne 3 ]] && [[ $(awk 'BEGIN{FS="\t"}{print NF}' $DATA | uniq) -ne 6 ]] && [[ $(awk 'BEGIN{FS="\t"}{print NF}' $DATA | uniq) -ne 12 ]]; then
		printf "%s\n" "" "ERROR: The structure of your file (${DATA}) does not present the classical structure of bed file" "This script has been developed to work this bed3, bed6 or bed12" "Please try with an orther file"
		exit 1
	else
		awk 'BEGIN{FS="\t"; OFS="\t"}{print ($3-$2)+1}' $DATA | awk 'BEGIN{FS="\t"; OFS="\t"}{SUM += $1} END{print (SUM/"'$BASEHUMANGENOME'")*100}' > .pourcent.tmp & PID=$!
		i=0 &
		while kill -0 $PID 2>/dev/null; do
			i=$(( (i+1) %4 ))
			printf "\rCalculation in process ${SPIN:$i:1}"
		sleep .1
		done
		POURCENTFILE="$(cat .pourcent.tmp)"
		rm .pourcent.tmp
		printf "\n\nThe coverage of the regions in your BED file represents "$POURCENTFILE"%% of the human genome GRCh37p13 (hg19).\n\n"
		exit 0
	fi	
elif [[ $(awk 'BEGIN{FS="\t"}{print NF}' $DATA | uniq) -eq 3 ]] || [[ $(awk 'BEGIN{FS="\t"}{print NF}' $DATA | uniq) -eq 6 ]] || [[ $(awk 'BEGIN{FS="\t"}{print NF}' $DATA | uniq) -eq 12 ]]; then
	awk 'BEGIN{FS="\t"; OFS="\t"}{print ($3-$2)+1}' $DATA | awk 'BEGIN{FS="\t"; OFS="\t"}{SUM += $1} END{print (SUM/"'$BASEHUMANGENOME'")*100}' > .pourcent.tmp & PID=$!
	i=0 &
	while kill -0 $PID 2>/dev/null; do
		i=$(( (i+1) %4 ))
		printf "\rCalculation in process ${SPIN:$i:1}"
	sleep .1
	done
	POURCENTFILE="$(cat .pourcent.tmp)"
	rm .pourcent.tmp
	printf "\n\nThe coverage of the regions in your BED file represents "$POURCENTFILE"%% of the human genome GRCh37p13 (hg19).\n\n"
	exit 0
fi
printf "%s\n" "" "Your file (${DATA}) does not present the classical structure of .gff or .bed file" "Please try with an orther file"
exit 1