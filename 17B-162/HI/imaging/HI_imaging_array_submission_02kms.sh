#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --mem=128000M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --job-name=M33_bc_02kms-%A-%a
#SBATCH --output=casa-m33_bc_02kms-%A-%a.out
#SBATCH --array=0-334

# Use array to set which channels will get imaged.
# Run from a separate folder so the log files are in one place.

module restore my_default

source /home/ekoch/.bashrc
source /home/ekoch/preload.bash

job_num=$SLURM_ARRAY_TASK_ID

# Build in a slight offset for each job to avoid starting a bunch of CASA
# sessions at once
sleep $(($job_num + 30))

# Move to scratch space b/c casa write out the temporary files into the same folder
export scratch_path=/home/ekoch/scratch/17B-162_imaging/

cd $scratch_path

Xvfb :1 &
export DISPLAY=:1

# Set which imaging stage to run: the initial 5-sigma clean (1) or the full clean with a mask to (2)
stage=1

if (( $stage==1 )); then
    script_name="${HOME}/code/VLA_Lband/17B-162/HI/imaging/HI_single_channel_clean.py"
    # Parameter file for tclean
    param_file="/home/ekoch/code/VLA_Lband/17B-162/HI/imaging/param_files/14B_17B_02kms.saved"
elif (( $stage==2 )); then
    script_name="${HOME}/code/VLA_Lband/17B-162/HI/imaging/HI_single_channel_clean_stage2.py"
    param_file="/home/ekoch/code/VLA_Lband/17B-162/HI/imaging/param_files/14B_17B_02kms_stage2.saved"
else
    echo "Stage must be 1 or 2, not ${stage}".
    exit 1
fi

# Start 5 channels running on the node
# This is well-suited for the cedar base nodes
start_chan=$(($job_num * 4))
end_chan=$((($job_num + 1) * 4))

# B/c CASA spawns other processes internally, and if something crashes in the scripts,
# the wait command may hang until the job is killed
# Try recording the casa interpreter PIDs and only make wait subject to those.
pids=

for (( chan_num = $start_chan; chan_num < $end_chan; chan_num++ )); do

    echo "Running channel "$chan_num

    $HOME/casa-release-5.3.0-143.el7/bin/mpicasa -n 8 $HOME/casa-release-5.3.0-143.el7/bin/casa --nologger --nogui --logfile $scratch_path/HI_contsub_021kms/casa_M33_HI_14B_17B_02kms_${chan_num}_${SLURM_JOB_ID}_$(date "+%Y%m%d-%H%M%S")_stage${stage}.log --nocrashreport -c $script_name $chan_num $param_file "HI_contsub_021kms" &
    pids+=" $!"

    # Avoid starting up multiple CASA sessions at once. Runs into ipython database issues.
    sleep 60

done

wait $pids || { echo "There was an error" >&2; exit 1; }

echo "All CASA jobs exited."
