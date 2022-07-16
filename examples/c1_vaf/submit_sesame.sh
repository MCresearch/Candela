#!/bin/bash
#PBS -l nodes=1:ppn=1,walltime=24:00:00
###PBS -q risk
#PBS -q eac 
#PBS -r n
#PBS -N VELCOR

cd $PBS_O_WORKDIR
/home/mohanc/3_DIY_Programs/D310/src/D310.exe > Log.txt
exit $?

mkdir DONE
