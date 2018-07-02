#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=128000M
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH --job-name=M33_dirty_cube-%J
#SBATCH --output=casa-m33_dirtycube-%J.out
export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

module restore my_default

source /home/ekoch/.bashrc
source /home/ekoch/preload.bash

# Give which spw using the --export flag when calling sbatch

export scratch_path=/home/ekoch/scratch/17B-162_imaging/
# export project_path=/home/ekoch/project/ekoch/

# Move to scratch space b/c casa write out the temporary files into the same folder
cd $scratch_path

Xvfb :1 &
export DISPLAY=:1

echo "OMP_NUM_THREADS "$OMP_NUM_THREADS
echo "Running SPW "$spw

spw=5

$HOME/casa-release-5.3.0-143.el7/bin/mpicasa -n 64 $HOME/casa-release-5.3.0-143.el7/bin/casa --nologger --nogui --log2term --nocrashreport -c $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.py $spw

# Copy the dirty_cube folder into project space
# cp -R $scratch_path/dirty_cube $project_path
