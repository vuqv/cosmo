#!/bin/bash
## This bah script to run on default allocation
#PBS -N cat_SETINDEX
##PBS -r n
#PBS -o output.out
#PBS -e error.err
#PBS -l nodes=1:ppn=1
#PBS -l walltime=8:00:00:00
##PBS -j oe
#PBS -A cyberlamp -l qos=cl_open
cd $PBS_O_WORKDIR

cur_dir=`pwd`
source /storage/home/qvv5013/work/anaconda3/etc/profile.d/conda.sh
#conda activate py37
python single_run_get_time_fQ.py -f control_cal_fQ.config