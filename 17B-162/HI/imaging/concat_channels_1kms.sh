#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=4000M
#SBATCH --ntasks=1
#SBATCH --job-name=M33_HI_concat_channel-1kms-%A-%a
#SBATCH --output=casa-m33_HI_concat_channel-1kms-%A-%a.out
#SBATCH --array=0-9

# Combine the imaged channels into cubes
# Run on whole node for the memory.

module restore my_default

source /home/ekoch/.bashrc
source /home/ekoch/preload.bash

export scratch_path=/home/ekoch/scratch/17B-162_imaging/

# Move to scratch space b/c casa write out the temporary files into the same folder
# cd $scratch_path

Xvfb :1 &
export DISPLAY=:1

suffixes="mask model pb psf residual residual_init image sumwt weight"
suffix_arr=($suffixes)

job_num=$SLURM_ARRAY_TASK_ID

# CASA runs better if each job has it's own source files
# So we'll run from the node storage
tmp_dir=$SLURM_TMPDIR/concat_M33_chans_1kms_${suffix_arr[$job_num]}
mkdir $tmp_dir

cd $tmp_dir

export casa_scratch_path="$HOME/scratch/casa-release-5.4.1-32.el7"

cp -r $casa_scratch_path .

# Copy the init file
mkdir .casa
cp $HOME/.casa/init.py .casa/

rc_path="${tmp_dir}/.casa"

# Now copy all of the images to the node storage
mkdir HI_contsub_1_0kms
for i in {0..272}; do
    mkdir HI_contsub_1_0kms/channel_${i}
    cp -r $scratch_path/HI_contsub_1_0kms/channel_${i}/*.${suffix_arr[$job_num]} HI_contsub_1_0kms/channel_${i}/
    echo "Copied channel ${i}"
done

echo "Concatenating channels"
casa-release-5.4.1-32.el7/bin/casa --rcdir ${rc_path} --nologger --nogui --log2term --nocrashreport -c $HOME/code/VLA_Lband/17B-162/HI/imaging/concat_channels.py "HI_contsub_1_0kms" "M33_14B_17B_HI_contsub_width_1kms" 272 ${suffix_arr[$job_num]}

echo "Copying M33_14B_17B_HI_contsub_width_1kms.${suffix_arr[$job_num]} to scratch"
cp -r "HI_contsub_1_0kms/M33_14B_17B_HI_contsub_width_1kms.${suffix_arr[$job_num]}" $scratch_path/HI_contsub_1_0kms

echo "Done!"