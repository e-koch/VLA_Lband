#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH --mem=256GB
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --job-name=M33_turbulence
#SBATCH --output=m33_turbulence-%J.out
export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

module restore my_default

# Use version loaded while building pyfftw
module load fftw-mpi/3.3.6

source /home/ekoch/.bashrc
source /home/ekoch/preload.bash

export scratch_path=/home/ekoch/scratch/14B-088
export project_path=/home/ekoch/project/ekoch/

# Call script with number of cores
$HOME/anaconda3/bin/python $HOME/code/VLA_Lband/14B-088/turbulence/M33_turbstats.py 32

# Copy the output files to the project path
cp $scratch_path/M33_turbulence/*.pkl $project_path/M33_turbulence/
