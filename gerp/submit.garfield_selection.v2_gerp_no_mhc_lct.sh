#!/bin/bash
  
#PBS -N garfield-gerp   
#PBS -o garfield-gerp-output
#PBS -e garfield-gerp-error
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=8
#PBS -S /bin/bash
# PBS -t 1-22

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

set -e

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

i=${PBS_ARRAYID}

cd /panfs/panasas01/shared-godmc/GARFIELDv2/garfield-v2

echo gerp
./garfield_gerp_no_mhc_lct
echo gerp_cis
./garfield_gerp_cis_no_mhc_lct
echo gerp_trans
./garfield_gerp_trans_no_mhc_lct
echo gerp_ambivalent
./garfield_gerp_ambivalent_no_mhc_lct

