#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=515000M
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH --job-name=M33_16B-242_dirty_cube-%A-%a
#SBATCH --output=casa-m33_16B-242_dirtycube-%A-%a.out
#SBATCH --array=0-9%1

export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

module restore my_default

source /home/ekoch/.bashrc
source /home/ekoch/preload.bash

# Give which spw using the --export flag when calling sbatch

export scratch_path=/home/ekoch/scratch/16B-242_imaging/

# Move to scratch space b/c casa write out the temporary files into the same folder
cd $scratch_path

Xvfb :1 &
export DISPLAY=:1

echo "OMP_NUM_THREADS "$OMP_NUM_THREADS
echo "Running SPW "$spw

spw=$SLURM_ARRAY_TASK_ID

$HOME/casa-release-5.3.0-143.el7/bin/mpicasa -n 32 $HOME/casa-release-5.3.0-143.el7/bin/casa --nologger --nogui --log2term --nocrashreport -c $HOME/code/VLA_Lband/16B/16B-242/imaging/test_line_imaging.py $spw
