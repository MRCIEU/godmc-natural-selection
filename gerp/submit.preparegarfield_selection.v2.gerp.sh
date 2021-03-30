#!/bin/bash
  
#PBS -N prepare-garfield   
#PBS -o garfieldprep-output
#PBS -e garfieldprep-error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=4
#PBS -S /bin/bash
#PBS -t 1-22

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

set -e

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

cd ~/repo/godmc_phase2_analysis/13_selection

R CMD BATCH --no-save --no-restore '--args '$i'' preparegarfield_selection_no_mhc_lct_gerp.R preparegarfield_selection_no_mhc_lctgerp$i.Rout



#BEFORE RUNNING THIS SCRIPT:
# you need to run 
# preparegarfield_selection0.R
# compile.snpcpgpval.sh
