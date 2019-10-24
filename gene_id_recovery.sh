#!/bin/bash
#
#Made by Mathias Boulanger
#Version 1
#
#Use to extract the geneID of offical_named gene list from gff_file
#usage: gene_id_recovery.sh my_list my_gff

ARGS=2
NAMEPROG=$(basename ${0})
LIST=$1
NAMELIST=$(basename ${LIST})
GFF=$2
NAMEGFF=$(basename ${GFF})
EXTENSION="list"
SPIN='-\|/'
RED='\033[1;31m'
GREEN='\033[1;32m'
ORANGE='\033[0;33m'
NOCOLOR='\033[0m'

NAMEFILE_nomatch=${NAMELIST%%.*}_no_geneID_found.${EXTENSION}
NAMEFILE_match=${NAMELIST%%.*}_geneID.${EXTENSION}
NAMEFILE_match_geneIDonly=${NAMELIST%%.*}_geneID_only.${EXTENSION}
NAMEFILE_multiplematch=${NAMELIST%%.*}_multiple_geneID.${EXTENSION}

#checking files
if [[ $# -ne $ARGS ]]; then
    printf "\n${GREEN}Usage:${NOCOLOR} ${NAMEPROG} my_list my_gff\n\n"
    exit 1
elif [[ ! -f $LIST ]];then
	printf "\n${RED}Error:${NOCOLOR} the file '${LIST}' does not exit!\nPlease use an existing file.\n${GREEN}Usage:${NOCOLOR} ${NAMEPROG} my_list my_gff\n\n"
	exit 1
elif [[ $(wc -l $LIST) = "0 ${LIST}" ]]; then
	printf "\n${RED}Error:${NOCOLOR} the file '${LIST}' is empty!\n\n"
	exit 1
elif [[ ! -f $GFF ]];then
	printf "\n${RED}Error:${NOCOLOR} the file '${GFF}' does not exit!\nPlease use an existing file.\n${GREEN}Usage:${NOCOLOR} ${NAMEPROG} my_list my_gff\n\n"
	exit 1
elif [[ $(wc -l $GFF) = "0 ${GFF}" ]]; then
	printf "\n${RED}Error:${NOCOLOR} the file '${GFF}' is empty!\n\n"
	exit 1
fi
##Trash all tmp file if it is exist
rm -f /tmp/${NAMEPROG}_*.tmp

question () {
	while true; do
		printf "%s\n" "The file ${1} exists!" "Do you want to overwrite this file? (Y/n)" 
		read ANSWER
		printf "\n"
		case $ANSWER in
			[yY]|[yY][eE][sS]|"" )
				rm -f ${1}
				break;;
			[nN]|[nN][oO] )
				printf "%s\n" "script stopped by the user!" 
				exit 0
				break;;
			* )
				printf "\033c"
				printf "%s\n" "" "Please answer yes or no." ""
				;;
		esac
	done
}

if [[ -f $NAMEFILE_nomatch ]]; then
	question $NAMEFILE_nomatch
fi
if [[ -f $NAMEFILE_match ]]; then
	question $NAMEFILE_match
fi
if [[ -f NAMEFILE_match_geneIDonly ]]; then
	question $NAMEFILE_multiplematch
fi
if [[ -f $NAMEFILE_multiplematch ]]; then
	question $NAMEFILE_multiplematch
fi

#recovery of gff attributes
printf "\n\n"
cut -f9 $GFF | sed -e 's/\;/	/g' > /tmp/${NAMEPROG}_attributes0.tmp & PID=$!
i=0 &
while kill -0 $PID 2>/dev/null; do
	i=$(( (i+1) %4 ))
	printf "\rExtracting attributes present in the gff file ${SPIN:$i:1}"
	sleep .1
done
printf "\n\n"

#extraction of the gene list
NUMGENE=$(cat ${LIST} | wc -l | sed 's/ //g')
if [[ $NUMGENE -eq 1 ]]; then
	printf "The list contains only one gene!\n\n"
else
	printf "The list contains ${NUMGENE} genes!\n\n"
fi

#recovery of the 
for (( i = 1; i < $NUMGENE+1; i++ )); do
	eval gene_name="$(sed -n $i'p' ${LIST})"
	cat /tmp/${NAMEPROG}_attributes0.tmp | grep "gene=${gene_name}\t" > /tmp/${NAMEPROG}_attributes1.tmp
	NUM_MATCH_ATTRIBUTE=$(cat /tmp/${NAMEPROG}_attributes1.tmp | wc -l)
	if [[ $NUM_MATCH_ATTRIBUTE -eq 0 ]]; then
		printf "${gene_name}\n" >> $NAMEFILE_nomatch
	elif [[ $NUM_MATCH_ATTRIBUTE -ge 2  ]]; then		
		if [[ $(awk 'BEGIN{FS="\t"}{n[$2]++}END{for (x in n) {print n[x]}}' /tmp/${NAMEPROG}_attributes1.tmp | sed -n '1p') -eq $NUM_MATCH_ATTRIBUTE ]]; then
			while true; do
				for (( k = 1; k < $(sed -n '1p' /tmp/${NAMEPROG}_attributes1.tmp | awk 'BEGIN{FS="\t"}{print NF}') +1; k++ )); do
					sed -n '1p' /tmp/${NAMEPROG}_attributes1.tmp | awk 'BEGIN{FS="\t";OFS="\t"}{split($'$k', subfield, "="); if(subfield[1]=="Dbxref") print "'${gene_name}'", subfield[2]}' >> $NAMEFILE_match
					if [[ $(cat $NAMEFILE_match | wc -l) -eq $i ]]; then
						break 2
					fi
				done
			done
		else
			printf "${gene_name}\n" >> $NAMEFILE_multiplematch
		fi		
	else
		while true; do
			for (( k = 1; k < $(awk 'BEGIN{FS="\t"}{print NF}' /tmp/${NAMEPROG}_attributes1.tmp) +1; k++ )); do
				awk 'BEGIN{FS="\t";OFS="\t"}{split($'$k', subfield, "="); if(subfield[1]=="Dbxref") print "'${gene_name}'", subfield[2]}' /tmp/${NAMEPROG}_attributes1.tmp >> $NAMEFILE_match
				if [[ $(cat $NAMEFILE_match | wc -l) -eq $i ]]; then
					break 2
				fi
			done
		done
	fi
	printf "\rRecovery of the gene num ${i} on ${NUMGENE} "
done

printf "\n"
if [[ -f $NAMEFILE_nomatch ]]; then
	NUM_NOMATCH=$(cat $NAMEFILE_nomatch | wc -l)
	if [[ $NUM_NOMATCH -eq 1 ]]; then
		printf "\n${ORANGE}Warning:${NOCOLOR} 1 gene of the list has not been found in the gff file!\n\n"
	else
		printf "\n${ORANGE}Warning:${NOCOLOR} ${NUM_NOMATCH} genes of the list have not been found in the gff file!\n\n"
	fi
fi
if [[ -f $NAMEFILE_multiplematch ]]; then
	NUM_MULTIPLE_MATCH=$(cat $NAMEFILE_multiplematch | wc -l)
	if [[ $NUM_MULTIPLE_MATCH -eq 1 ]]; then
		printf "\n${ORANGE}Warning:${NOCOLOR} 1 gene of the list has been found multiple time with different geneIDs in the gff file!\n\n"
	else
		printf "\n${ORANGE}Warning:${NOCOLOR} ${NUM_MULTIPLE_MATCH} gene of the list have been found multiple time with different geneIDs in the gff file!\n\n"
	fi	
fi

#creation of ID list
if [[ -f $NAMEFILE_match ]]; then
	cut -f2 $NAMEFILE_match | awk 'BEGIN{FS=":";OFS="\t"}{split($2, subfield1, ","); print subfield1[1]}' > $NAMEFILE_match_geneIDonly
fi

printf "%s\n" "" "GeneIDs have been extracted with success!" ""
rm /tmp/${NAMEPROG}_attributes0.tmp
rm /tmp/${NAMEPROG}_attributes1.tmp

exit 0
