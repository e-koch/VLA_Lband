#!/bin/bash
#SBATCH -c 1
#SBATCH --time=12:00:00
#SBATCH --mem=16000M
#SBATCH --job-name=M33_dirty_cube-%A-%a
#SBATCH --output=casa-m33_dirtycube-%A-%a.out
#SBATCH --array=2400-2410

# Use array to set which channels will get imaged.
# Run from a separate folder so the log files are in one place.

module restore my_default

source /home/ekoch/.bashrc
source /home/ekoch/preload.bash

chan_num=$SLURM_ARRAY_TASK_ID

# Move to scratch space b/c casa write out the temporary files into the same folder
export scratch_path=/home/ekoch/scratch/17B-162_imaging/

cd $scratch_path

Xvfb :1 &
export DISPLAY=:1

echo "Running channel "$chan_num

$HOME/casa-release-5.3.0-143.el7/bin/casa --nologger --nogui --log2term --nocrashreport -c $HOME/code/VLA_Lband/17B-162/HI/imaging/HI_single_channel_clean.py $chan_num

