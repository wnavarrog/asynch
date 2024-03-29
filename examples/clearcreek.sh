#!/bin/sh
#$ -N clearcreek_run
#$ -j y
#$ -cwd
####$ -m e
####$ -M your_email_address_goes_here
#$ -pe orte 8
####$ -pe 8cpn 8
####$ -l mf=16G
#$ -l ib=1
####$ -q UI
#$ -q IFC
####$ -q all.q
####$ -q COE

/bin/echo Running on host: `hostname`.
/bin/echo In directory: `pwd`
/bin/echo Starting on: `date`
/bin/echo "Got $NSLOTS processors."
/bin/echo 
/bin/echo 
/bin/echo 

mpirun -np 8 ./ASYNCH examples/Global254.gbl

