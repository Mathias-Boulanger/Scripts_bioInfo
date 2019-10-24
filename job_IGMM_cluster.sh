#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -N Job_name
#$ -V
#$ -o Log_name.log
#$ -M mathias.boulanger@igmm.cnrs.fr
#$ -S /bin/bash
#
#Launch this script like this: qsub -q {qshort|qlong} -o log.txt job_cluster.sh

printf "%s\n" "" "----------------------------------------" "-----------Starting job date:-----------"
date
printf "%s\n" "" "----------------------------------------" "--------------Launched by:--------------"
whoami
printf "%s\n" "" "----------------------------------------" "-----------------JOB ID:----------------" "${JOB_ID}"
printf "%s\n" "" "----------------------------------------" "-------------execution Node:------------" "${HOSTNAME}"
printf "%s\n" "" "----------------------------------------" ""

#export in the path the anaconda3 MPL environnement
export PATH="/homegalio2/mpl/Tools/anaconda3/bin:$PATH"

#Check which tool and which version
which TOOLS
TOOLS --version

#Job command:




#Job end
printf "%s\n" "" "----------------------------------------" "------------Ending job date:------------"
date
printf "%s\n" "" "----------------------------------------" ""
