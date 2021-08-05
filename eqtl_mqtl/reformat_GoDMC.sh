#!/bin/bash

#PBS -N godmc
#PBS -o godmc-output
#PBS -e godmc-error
# PBS -t 1-962
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1
#PBS -q himem
#PBS -S /bin/bash

cd /newhome/epzjlm/repo/godmc-natural-selection/eqtl_mqtl/
Rscript reformat_GoDMC.R
