#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=515000M
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH --job-name=M33_16B-236_clean_cube-%J
#SBATCH --output=casa-m33_16B-236_cleancube-%J.out

export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

module restore my_default

source /home/ekoch/.bashrc
source /home/ekoch/preload.bash

# Give which spw using the --export flag when calling sbatch

export scratch_path=/home/ekoch/scratch/16B-236_imaging/

# Move to scratch space b/c casa write out the temporary files into the same folder
cd $scratch_path

Xvfb :1 &
export DISPLAY=:1

echo "OMP_NUM_THREADS "$OMP_NUM_THREADS

$HOME/casa-release-5.3.0-143.el7/bin/mpicasa -n $OMP_NUM_THREADS $HOME/casa-release-5.3.0-143.el7/bin/casa --nologger --nogui --log2term --nocrashreport -c $HOME/code/VLA_Lband/16B/16B-236/HI/imaging/HI_cube_imaging.py
