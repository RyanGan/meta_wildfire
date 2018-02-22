#!/bin/bash
#PBS -N ryan_rjob
#PBS -l nodes=1:ppn=6
#PBS -W group_list=pierce_group
#PBS -q batch
#PBS -M rgan@colostate.edu
#PBS -m abe

# set memory to unlimited
ulimit -s unlimited

# set command directory
cd /home/ryangan/local_git_repo/meta_wildfire/

# export R library
export R_LIB$=$R_LIB:/home/ryangan/R/x86_64-pc-linux-gnu-library/3.4/

# run rscript
Rscript --vanilla ./r_scripts/exposure_wrangle/1016-narr_temp

