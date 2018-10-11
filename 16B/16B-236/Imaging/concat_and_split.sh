#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --mem=128000M
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH --job-name=M33_236_concat_split-%J
#SBATCH --output=casa-m33_236_concat_split-%J.out

export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

module restore my_default

source /home/ekoch/.bashrc
source /home/ekoch/preload.bash

export scratch_path=/home/ekoch/scratch/16B-236_imaging/

# Move to scratch space b/c casa writes out the temporary files into the same folder
cd $scratch_path

Xvfb :1 &
export DISPLAY=:1

$HOME/casa-release-5.3.0-143.el7/bin/casa --nologger --nogui --log2term --nocrashreport -c $HOME/code/VLA_Lband/16B/16B-236/Imaging/concat_and_split.py
