#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
#Launch this script like this: qsub -q allnodes -o log.txt my_script.sh

#export user PATH to the cluster 
#export PATH="$SGE_O_PATH"

#print d

bowtie2-build /homegalio2/boulangerm/GRCh37p13_hg19_MB/GRCh37_20190212_mainChr_chrFormat.fna GRCh37p13_20190212_mainChr_chrFormat

exit 0