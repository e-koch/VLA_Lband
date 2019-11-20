#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=128000M
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH --job-name=M33_HI_match_and_combine-%A-%a
#SBATCH --output=casa-m33_HI_match_and_combine-%A-%a.out
#SBATCH --array=0-9%1


# Create regridded 14B and 17B HI MMS and split into individual channels

export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

module restore my_default

source /home/ekoch/.bashrc
source /home/ekoch/preload.bash

export scratch_path=/home/ekoch/scratch/17B-162_imaging/

job_num=$SLURM_ARRAY_TASK_ID

# CASA runs better if each job has it's own source files
# So we'll run from the node storage
tmp_dir=$SLURM_TMPDIR/concat_M33_chans_04kms_${suffix_arr[$job_num]}
mkdir $tmp_dir

cd $tmp_dir

export casa_path="$HOME/casa-release-5.4.1-32.el7"

cp -r $casa_path .

# Copy the init file
mkdir .casa
cp $HOME/.casa/init.py .casa/

rc_path="${tmp_dir}/.casa"

Xvfb :1 &
export DISPLAY=:1

export out_chan_folder="HI_nocontsub_0_42kms"

# Move to scratch space b/c casa write out the temporary files into the same folder
# cd $scratch_path

echo "Copy whole MS to node"

# cp -r $scratch_path/17B-162_HI_spw_0_LSRK.ms.contsub .
# cp -r $scratch_path/14B-088_HI_LSRK.ms.contsub .

# Run match_and_split_chunks.py to create smaller chunked versions of the MSs first.
# The 17B data is too large to fit on node storage.

# Only copy chunk MMSs over.
cp -r $scratch_path/17B-162_HI_spw_0_LSRK.ms.0_42kms.mms_chunk_${job_num} .
cp -r $scratch_path/14B-088_HI_LSRK.ms.0_42kms.mms_chunk_${job_num} .

mkdir ${out_chan_folder}

casa-release-5.4.1-32.el7/bin/mpicasa -n 32 casa-release-5.4.1-32.el7/bin/casa --rcdir ${rc_path} --nologger --nogui --log2term --nocrashreport -c $HOME/code/VLA_Lband/17B-162/HI/imaging/match_and_split.py False 0.42 $job_num 10 ${scratch_path}

# Local run for 0.42 km/s. NOTE THE DIFFERENCE
# ~/casa-release-5.4.1-32.el7/bin/casa --nologger --nogui --log2term --nocrashreport -c $HOME/ownCloud/code_development/VLA_Lband/17B-162/HI/imaging/match_and_split.py False 0.42
