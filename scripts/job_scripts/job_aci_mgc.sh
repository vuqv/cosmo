#!/bin/bash
## This bah script to run on mgc-nih allocation
#PBS -N N200_CONFIDX
#PBS -l walltime=24:00:00
#PBS -A mgc_nih
#PBS -l qos=mgc_nih
#PBS -l pmem=2gb
#PBS -j oe
#PBS -l feature=rhel7
#PBS -l nodes=1:ppn=1:gpus=1:shared:gc_t4

cd $PBS_O_WORKDIR
echo `pwd`

source /storage/home/qvv5013/work/anaconda3/etc/profile.d/conda.sh

python single_run_get_time_fQ.py -f control_cal_fQ.config

