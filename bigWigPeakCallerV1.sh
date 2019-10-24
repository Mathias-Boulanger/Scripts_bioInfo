#!/bin/bash
#
#Made by Mathias Boulanger and Fabienne Bejjani
#Version BETA
#
#Need biWigToWig and WigPeakCaller utilities in the same folder than the present script.
#USAGE: bash bigWigPeakCaller.sh bigWigFile.bw threshold minRun maxGap

ARGS=4						#The script need X argument
NAMEPROG=$(basename ${0})	#Name of the program
pathToProg=$(dirname ${0})
BIGWIG=$1					#File in argument
NAMEBIGWIG=$(basename ${BIGWIG})
pathToBigWig=$(dirname ${BIGWIG})
WORKDIR=${pathToBigWig}/${NAMEBIGWIG%%.*}_${NAMEPROG%%.*}
EXTENSION1="bw"				#File extension 1  necessary to run this script
EXTENSION2="bigWig"			#File extension 2  necessary to run this script
SPIN='-\|/'					#Waiting characters
RED='\033[1;31m'
GREEN='\033[1;32m'
ORANGE='\033[0;33m'
NOCOLOR='\033[0m'

##Check the ability to work
#Check the presence and the executability of the required elements
if [[ ! -f ${pathToProg}/bigWigToWig ]]; then
	printf "\n${RED}Error:${NOCOLOR} The bigWigToWig utility is not in the same directory than ${NAMEPROG}!\n\n"
	exit 1
elif [[ ! -x ${pathToProg}/bigWigToWig ]]; then
	printf "\n${RED}Error:${NOCOLOR} The bigWigToWig utility is not executable!\n\n"
	exit 1
elif [[ ! -f ${pathToProg}/WigPeakCaller ]]; then
	printf "\n${RED}Error:${NOCOLOR} The WigPeakCaller utility is not in the same directory than ${NAMEPROG}!\n\n"
	exit 1
elif [[ ! -x ${pathToProg}/WigPeakCaller ]]; then
	printf "\n${RED}Error:${NOCOLOR} The WigPeakCaller utility is not executable!\n\n"
	exit 1
fi

#check input arguments
if [[ $# -ne $ARGS ]]; then
	printf "\n${RED}Error:${NOCOLOR} Not correct number of argument (${#})\n${GREEN}Usage:${NOCOLOR} bash ${NAMEPROG} bigWigFile.bw threshold minRun maxGap\n\n"
	exit 1
elif [[ ! -f $BIGWIG ]];then
	printf "\n${RED}Error:${NOCOLOR} the file '${BIGWIG}' does not exit!\nPlease use an existing file.\n${GREEN}Usage:${NOCOLOR} bash ${NAMEPROG} bigWigFile.bw threshold minRun maxGap\n\n"
	exit 1
elif [[ $(wc -l $BIGWIG) = "0 ${BIGWIG}" ]]; then
	printf "\n${RED}Error:${NOCOLOR} the file is empty!\n\n"
	exit 1
elif [[ ! "${2}" =~ ^[-+]?[0-9]+$ ]]; then
	printf "\n${RED}Error:${NOCOLOR} The argument for the threshold is not numeric!\n${GREEN}Usage:${NOCOLOR} bash ${NAMEPROG} bigWigFile.bw threshold minRun maxGap\n\n"
	exit 1
elif [[ ! "${3}" =~ ^[-+]?[0-9]+$ ]]; then
	printf "\n${RED}Error:${NOCOLOR} The argument for the minRun is not numeric!\n${GREEN}Usage:${NOCOLOR} bash ${NAMEPROG} bigWigFile.bw threshold minRun maxGap\n\n"
	exit 1
elif [[ ! "${4}" =~ ^[-+]?[0-9]+$ ]]; then
	printf "\n${RED}Error:${NOCOLOR} The argument for the maxGap is not numeric!\n${GREEN}Usage:${NOCOLOR} bash ${NAMEPROG} bigWigFile.bw threshold minRun maxGap\n\n"
	exit 1
fi




#Check the extention of the inupt file
if [[ ${BIGWIG##*.} != $EXTENSION1 ]]; then
	if [[ ${BIGWIG##*.} != $EXTENSION2 ]]; then
		printf "\n${ORANGE}WARNING:${NOCOLOR} The file extension should be .${EXTENSION1} or .${EXTENSION2}\nMake sure that the file present is a bigwig file...\n"
	fi
fi

#Trash all tmp file if it is exist
rm -f /tmp/${NAMEPROG}_*.tmp

printf "\nReady to work!\n"
##Script Work
#work directory
if [[ ! -d ${WORKDIR} ]]; then
	mkdir ${WORKDIR}
fi

#Transform bigWig to wig
${pathToProg}/bigWigToWig ${BIGWIG} ${WORKDIR}/${NAMEBIGWIG%%.*}.wig & PID=$!
i=0 &
while kill -0 $PID 2>/dev/null; do
	i=$(( (i+1) %4 ))
	printf "\rbigWig to wig conversion  ${SPIN:$i:1}"
	sleep .1
done
printf "\n"

#Peak calling
${pathToProg}/WigPeakCaller --wigFolder ${WORKDIR} --wigRegex ${NAMEBIGWIG%%.*}.wig --threshold $2 --minRun $3 --maxGap $4 --outputFolder ${WORKDIR} & PID=$!
i=0 &
while kill -0 $PID 2>/dev/null; do
	i=$(( (i+1) %4 ))
	printf "\rPeak calling  ${SPIN:$i:1}"
	sleep .1
done
printf "\n"

printf "\n${GREEN}Done!${NOCOLOR}\n"
rm -f /tmp/${NAMEPROG}_*.tmp
exit 0