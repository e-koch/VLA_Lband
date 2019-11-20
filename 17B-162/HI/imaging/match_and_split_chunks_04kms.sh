#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=128000M
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH --job-name=M33_HI_match_and_combine_chunks-%J
#SBATCH --output=casa-m33_HI_match_and_combine_chunks-%J.out


# Create chunked 14B and 17B HI MMS to make splitting individual channels easier + faster

export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

module restore my_default

source /home/ekoch/.bashrc
source /home/ekoch/preload.bash

export scratch_path=/home/ekoch/scratch/17B-162_imaging/

cd $scratch_path

Xvfb :1 &
export DISPLAY=:1

export out_chan_folder="HI_nocontsub_0_42kms"

~/casa-release-5.4.1-32.el7/bin/mpicasa -n 32 ~/casa-release-5.4.1-32.el7/bin/casa --nologger --nogui --log2term --nocrashreport -c $HOME/code/VLA_Lband/17B-162/HI/imaging/match_and_split_chunks.py False 0.42 10
