#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=128000M
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH --job-name=M33_HI_concat_channel-02kms-%J
#SBATCH --output=casa-m33_HI_concat_channel-02kms-%J.out

# Combine the imaged channels into cubes
# Run on whole node for the memory.

# export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

module restore my_default

source /home/ekoch/.bashrc
source /home/ekoch/preload.bash

export scratch_path=/home/ekoch/scratch/17B-162_imaging/

# Move to scratch space b/c casa write out the temporary files into the same folder
cd $scratch_path

Xvfb :1 &
export DISPLAY=:1

$HOME/casa-release-5.3.0-143.el7/bin/casa --nologger --nogui --log2term --nocrashreport -c $HOME/code/VLA_Lband/17B-162/HI/imaging/concat_channels.py "HI_contsub_021kms" "M33_14B_17B_HI_contsub_width_02kms" 1220
