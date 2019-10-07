#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=256000M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --job-name=M33_17B_OH-%A-%a
#SBATCH --output=casa-m33_17B_OH-%A-%a.out
#SBATCH --array=0-7

# Use array to set which channels will get imaged.
# Run from a separate folder so the log files are in one place.

module restore my_default

source /home/ekoch/.bashrc
source /home/ekoch/preload.bash

job_num=$SLURM_ARRAY_TASK_ID

# Build in a slight offset for each job to avoid starting a bunch of CASA
# sessions at once
# sleep $(($job_num + 30))

# Move to scratch space b/c casa write out the temporary files into the same folder
export scratch_path=/home/ekoch/scratch/17B-162_imaging/

cd $scratch_path

Xvfb :1 &
export DISPLAY=:1

# Path to the casa files
export casa_scratch_path="$HOME/casa-release-5.5.0-149.el7"

# Make a new directory on the node storage
tmp_dir=$SLURM_TMPDIR/OH_imaging_${chan_num}
mkdir $tmp_dir

cd $tmp_dir

# Copy a new casa instance to avoid slower i/o on scratch or in home
cp -r $casa_scratch_path .

# Copy the init file
mkdir .casa
cp $HOME/.casa/init.py .casa/

rc_path="${tmp_dir}/.casa"

weighting='natural'

if (( $job_num==0 )); then

    # OH1612 contsub
    mkdir $scratch_path/OH1612_imaging

    # Copy data
    cp -r $scratch_path/17B-162_OH1612_spw_3_LSRK.ms.contsub .

    casa-release-5.5.0-149.el7/bin/casa --rcdir ${rc_path} --nologger --nogui --log2term --logfile $scratch_path/OH_17B_imaging_$(date "+%Y%m%d-%H%M%S")_1612.log --nocrashreport -c ~/code/VLA_Lband/17B-162/OH/OH_imaging.py 3 True $weighting False

    cp -r OH1612_imaging/* $scratch_path/OH1612_imaging/

elif (( $job_num==1 )); then

    # OH1612 no contsub
    mkdir $scratch_path/OH1612_imaging

    # Copy data
    cp -r $scratch_path/17B-162_OH1612_spw_3_LSRK.ms .

    casa-release-5.5.0-149.el7/bin/casa --rcdir ${rc_path} --nologger --nogui --log2term --logfile $scratch_path/OH_17B_imaging_$(date "+%Y%m%d-%H%M%S")_1612_nocontsub.log --nocrashreport -c ~/code/VLA_Lband/17B-162/OH/OH_imaging.py 3 False $weighting False

    cp -r OH1612_imaging/* $scratch_path/OH1612_imaging/

elif (( $job_num==2 )); then

    # OH1665 contsub
    mkdir $scratch_path/OH1665_imaging

    # Copy data
    cp -r $scratch_path/17B-162_OH1665_spw_5_LSRK.ms.contsub .

    casa-release-5.5.0-149.el7/bin/casa --rcdir ${rc_path} --nologger --nogui --log2term --logfile $scratch_path/OH_17B_imaging_$(date "+%Y%m%d-%H%M%S")_1665.log --nocrashreport -c ~/code/VLA_Lband/17B-162/OH/OH_imaging.py 5 True $weighting False

    cp -r OH1665_imaging/* $scratch_path/OH1665_imaging/

elif (( $job_num==3 )); then

    # OH1665 no contsub
    mkdir $scratch_path/OH1665_imaging

    # Copy data
    cp -r $scratch_path/17B-162_OH1665_spw_5_LSRK.ms .

    casa-release-5.5.0-149.el7/bin/casa --rcdir ${rc_path} --nologger --nogui --log2term --logfile $scratch_path/OH_17B_imaging_$(date "+%Y%m%d-%H%M%S")_1665_nocontsub.log --nocrashreport -c ~/code/VLA_Lband/17B-162/OH/OH_imaging.py 5 False $weighting False

    cp -r OH1665_imaging/* $scratch_path/OH1665_imaging/

elif (( $job_num==4 )); then

    # OH1667 contsub
    mkdir $scratch_path/OH1667_imaging

    # Copy data
    cp -r $scratch_path/17B-162_OH1667_spw_6_LSRK.ms.contsub .

    casa-release-5.5.0-149.el7/bin/casa --rcdir ${rc_path} --nologger --nogui --log2term --logfile $scratch_path/OH_17B_imaging_$(date "+%Y%m%d-%H%M%S")_1667.log --nocrashreport -c ~/code/VLA_Lband/17B-162/OH/OH_imaging.py 6 True $weighting False

    cp -r OH1667_imaging/* $scratch_path/OH1667_imaging/

elif (( $job_num==5 )); then

    # OH1667 no contsub
    mkdir $scratch_path/OH1667_imaging

    # Copy data
    cp -r $scratch_path/17B-162_OH1667_spw_6_LSRK.ms .

    casa-release-5.5.0-149.el7/bin/casa --rcdir ${rc_path} --nologger --nogui --log2term --logfile $scratch_path/OH_17B_imaging_$(date "+%Y%m%d-%H%M%S")_1667_nocontsub.log --nocrashreport -c ~/code/VLA_Lband/17B-162/OH/OH_imaging.py 6 False $weighting False

    cp -r OH1667_imaging/* $scratch_path/OH1667_imaging/

elif (( $job_num==6 )); then

    # OH1712 contsub
    mkdir $scratch_path/OH1720_imaging

    # Copy data
    cp -r $scratch_path/17B-162_OH1720_spw_7_LSRK.ms.contsub .

    casa-release-5.5.0-149.el7/bin/casa --rcdir ${rc_path} --nologger --nogui --log2term --logfile $scratch_path/OH_17B_imaging_$(date "+%Y%m%d-%H%M%S")_1720.log --nocrashreport -c ~/code/VLA_Lband/17B-162/OH/OH_imaging.py 7 True $weighting False

    cp -r OH1720_imaging/* $scratch_path/OH1720_imaging/

elif (( $job_num==7 )); then

    # OH1720 no contsub
    mkdir $scratch_path/OH1720_imaging

    # Copy data
    cp -r $scratch_path/17B-162_OH1720_spw_7_LSRK.ms .

    casa-release-5.5.0-149.el7/bin/casa --rcdir ${rc_path} --nologger --nogui --log2term --logfile $scratch_path/OH_17B_imaging_$(date "+%Y%m%d-%H%M%S")_1720_nocontsub.log --nocrashreport -c ~/code/VLA_Lband/17B-162/OH/OH_imaging.py 7 False $weighting False

    cp -r OH1720_imaging/* $scratch_path/OH1720_imaging/

else
    echo "Stage must be 1 or 2, not ${stage}".
    exit 1
fi

echo "All CASA jobs exited."
