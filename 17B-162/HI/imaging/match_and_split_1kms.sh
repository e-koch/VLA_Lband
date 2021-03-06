#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=128000M
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH --job-name=M33_HI_match_and_combine-1kms-%A-%a
#SBATCH --output=casa-m33_HI_match_and_combine-1kms-%A-%a.out
#SBATCH --array=0-1%1

# In array mode, one job will run at a time. The match_and_split.py script will
# Check the stage that the previous job got to, then pick up from where the last one
# ended.

# Create a combined 14B and 17B HI MMS

export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

module restore my_default

source /home/ekoch/.bashrc
source /home/ekoch/preload.bash

export scratch_path=/home/ekoch/scratch/17B-162_imaging/

# Move to scratch space b/c casa write out the temporary files into the same folder
cd $scratch_path

Xvfb :1 &
export DISPLAY=:1

$HOME/casa-release-5.4.1-32.el7/bin/mpicasa -n 32 $HOME/casa-release-5.4.1-32.el7/bin/casa --nologger --nogui --log2term --nocrashreport -c $HOME/code/VLA_Lband/17B-162/HI/imaging/match_and_split.py True 1
