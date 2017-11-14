#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH --mem=128GB
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --job-name=M33_dirty_cube
#SBATCH --output=casa-m33_dirtycube-%J.out
export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

module restore my_default

source /home/ekoch/.bashrc
source /home/ekoch/preload.bash

export scratch_path=/home/ekoch/scratch/14B-088
export project_path=/home/ekoch/project/ekoch/

# Move to scratch space b/c casa write out the temporary files into the same folder
cd $scratch_path

$HOME/casa-release-5.1.1-5.el7/bin/mpicasa -n 32 $HOME/casa-release-5.1.1-5.el7/bin/casa --nologger --nogui --log2term --nocrashreport -c $HOME/code/VLA_Lband/14B-088/HI/imaging/cedar/HI_dirty_cube_synthimage.py $scratch_path

# Copy the dirty_cube folder into project space
cp -R $scratch_path/dirty_cube $project_path
