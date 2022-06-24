#!/bin/bash -l
#SBATCH -J cat_SETINDEX
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=5GB
#SBATCH --time=72:00:00 
#SBATCH -A plgribo3gpu
#SBATCH --gres=gpu
#SBATCH -p plgrid-gpu
#SBATCH --output=output.out
#SBATCH --error=error.err
cd $SLURM_SUBMIT_DIR
srun /bin/hostname

conda activate py310
python single_run_get_time_fQ.py -f control_cal_fQ.config

