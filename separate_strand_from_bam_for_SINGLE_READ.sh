#!/bin/bash
#
#Made by Mathias Boulanger
#Version 1
#
#separate strand of single read BAM aligment file into 2 different BAM files
#samtools is required for this script!
#Usage: bash separate_strand_from_bam_for_SINGLE_READ.sh BAM_file_to_split.bam

#0x0010 16 strand of the query (0 for forward; 1 for reverse strand)

ARGS=1
if [[ $# -ne $ARGS ]]; then
	printf "\n${RED}Error:${NOCOLOR} Not correct number of argument (${#})\n${GREEN}Usage:${NOCOLOR} bash ${NAMEPROG} BAM_file_to_split.bam\n\n"
	exit 1
fi

NAMEPROG=$(basename ${0})
name_BAM=$(basename ${1})
pathToBAM=$(cd $(dirname ${1}) && pwd -P)
EXTENSION="bam"	
SPIN='-\|/'
RED='\033[1;31m'
GREEN='\033[1;32m'
ORANGE='\033[0;33m'
NOCOLOR='\033[0m'


##Check the ability to work
#Check the samtools is installed
if [[ "$(command -v samtools)" == "" ]]; then
	printf "\n${RED}Error:${NOCOLOR} samtools command is missing to execute ${NAMEPROG} properly!\nPlease install it on your system to use ${NAMEPROG}\n\n"
	exit 1
else
	printf "VERSION:\n"
	samtools --version
	printf "\n"
fi

#check input argument
if [[ ! -f ${pathToBAM}/${name_BAM} ]];then
	printf "\n${RED}Error:${NOCOLOR} the file '${name_BAM}' does not exit!\nPlease use an existing file.\n${GREEN}Usage:${NOCOLOR} bash ${NAMEPROG} BAM_file_to_split.bam\n\n"
	exit 1
elif [[ $(wc -l ${pathToBAM}/${name_BAM}) = "0 ${pathToBAM}/${name_BAM}" ]]; then
	printf "\n${RED}Error:${NOCOLOR} the file '${name_BAM}' is empty!\n\n"
	exit 1
elif [[ ${name_BAM##*.} != ${EXTENSION} ]]; then
	printf "\n${ORANGE}WARNING:${NOCOLOR} The file extension should be .${EXTENSION}\nMake sure that the file present is a .${EXTENSION} file...\n"
fi

##Split BAM depending of the strand
#Forward
samtools view -F 0x10 -h ${pathToBAM}/${name_BAM} | samtools view -bS -> ${pathToBAM}/${name_BAM%%.*}_for.${EXTENSION} & PID=$!
i=0 &
while kill -0 $PID 2>/dev/null; do
	i=$(( (i+1) %4 ))
	printf "\rExtracting forward aligned reads ${SPIN:$i:1}"
	sleep .1
done
printf "\nDone!\n\n"

#Reverse
samtools view -f 0x10 -h ${pathToBAM}/${name_BAM} | samtools view -bS - > ${pathToBAM}/${name_BAM%%.*}_rev.${EXTENSION} & PID=$!
i=0 &
while kill -0 $PID 2>/dev/null; do
	i=$(( (i+1) %4 ))
	printf "\rExtracting reverse aligned reads ${SPIN:$i:1}"
	sleep .1
done
printf "\n${GREEN}Done!${NOCOLOR}\nBoth strand has been split into two different files!\n"
exit 0
