#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --mem=515000M
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH --job-name=M33_16B-242_transform_and_uvsub-%A-%a
#SBATCH --output=casa-m33_16B-242_transform_and_uvsub-%A-%a.out
#SBATCH --array=0-9%1

export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

module restore my_default

source /home/ekoch/.bashrc
source /home/ekoch/preload.bash

export scratch_path=/home/ekoch/scratch/16B-242_imaging/

# Move to scratch space b/c casa write out the temporary files into the same folder
cd $scratch_path

Xvfb :1 &
export DISPLAY=:1

spw=$SLURM_ARRAY_TASK_ID

$HOME/casa-release-5.4.1-32.el7/bin/mpicasa -n 32 $HOME/casa-release-5.4.1-32.el7/bin/casa --nologger --nogui --log2term --nocrashreport -c $HOME/code/VLA_Lband/16B/16B-242/imaging/transform_and_uvsub.py $spw
