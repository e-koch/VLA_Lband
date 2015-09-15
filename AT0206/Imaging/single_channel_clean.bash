#!/bin/bash

# Source bash profile
shopt -s expand_aliases
source /home/ekoch/.bash_profile

# Input
channel=${1}

mask_model_channel=$(($channel - 9))

echo $channel
echo $mask_model_channel

ms_name=M33_b_c_channel_${channel}.ms.tgz
mask_name=M33_newmask_channel_${mask_model_channel}.image.tgz
model_name=M33_model_channel_${mask_model_channel}.image.tgz

# Set username. Otherwise CASA crashes.
export USER='ekoch'

# Get certificate
getCert

mkdir -m 777  ${TMPDIR}/proc_${channel}

rm -rf /home/ekoch/canfar_scripts
git clone https://github.com/e-koch/canfar_scripts.git /home/ekoch/canfar_scripts

# echo "Mounting data"
# mountvofs --vospace vos:MWSynthesis/ --mountpoint ${TMPDIR}/vos --cache_dir ${TMPDIR}/vos_cache

echo "Copying files onto VM"

vcp vos:MWSynthesis/Arecibo/newmask_channels/${mask_name} ${TMPDIR}/proc_${channel}
echo "Done M33_mask"

vcp vos:MWSynthesis/Arecibo/model_channels/${model_name} ${TMPDIR}/proc_${channel}
echo "Done M33_model"

vcp vos:MWSynthesis/VLA/archival/single_channels/${ms_name} ${TMPDIR}/proc_${channel}
echo "Done MS Set"

# Unmount
# fusermount -u ${TMPDIR}/vos

cd ${TMPDIR}/proc_${channel}

# Unzip the model and mask
tar -zxf ${mask_name}
tar -zxf ${model_name}
tar -zxf ${ms_name}

# Delete zip files
rm ${mask_name}
rm ${model_name}
rm ${ms_name}

# Rename the files w/o the zipped ends
ms_name=${ms_name%.tgz}
mask_name=${mask_name%.tgz}
model_name=${model_name%.tgz}

echo "Running CASA"

echo "Show contents"
ls -al ${TMPDIR}/proc_${channel}

casapy --nogui -c /home/ekoch/canfar_scripts/img_pipe/archival_data/single_channel_clean.py ${ms_name} ${model_name} ${mask_name}

# Compress the clean output to upload to VOS
tar -zcf ${ms_name%.ms}.clean.image.tar.gz ${ms_name%.ms}.clean.image
tar -zcf ${ms_name%.ms}.clean.mask.tar.gz ${ms_name%.ms}.clean.mask
tar -zcf ${ms_name%.ms}.clean.model.tar.gz ${ms_name%.ms}.clean.model
tar -zcf ${ms_name%.ms}.clean.psf.tar.gz ${ms_name%.ms}.clean.psf
tar -zcf ${ms_name%.ms}.clean.residual.tar.gz ${ms_name%.ms}.clean.residual
tar -zcf ${ms_name%.ms}.clean.flux.tar.gz ${ms_name%.ms}.clean.flux
tar -zcf ${ms_name%.ms}.clean.flux.pbcoverage.tar.gz ${ms_name%.ms}.clean.flux.pbcoverage

# Now copy over the relevant infos

vmkdir vos:MWSynthesis/VLA/archival/test_bad_cleans/channel_${channel}

vcp ${TMPDIR}/proc_${channel}/casa*.log vos:MWSynthesis/VLA/archival/test_bad_cleans/channel_${channel}

vcp ${TMPDIR}/proc_${channel}/${ms_name%.ms}.clean.image.tar.gz vos:MWSynthesis/VLA/archival/test_bad_cleans/channel_${channel}/
vcp ${TMPDIR}/proc_${channel}/${ms_name%.ms}.clean.mask.tar.gz vos:MWSynthesis/VLA/archival/test_bad_cleans/channel_${channel}/
vcp ${TMPDIR}/proc_${channel}/${ms_name%.ms}.clean.model.tar.gz vos:MWSynthesis/VLA/archival/test_bad_cleans/channel_${channel}/
vcp ${TMPDIR}/proc_${channel}/${ms_name%.ms}.clean.psf.tar.gz vos:MWSynthesis/VLA/archival/test_bad_cleans/channel_${channel}/
vcp ${TMPDIR}/proc_${channel}/${ms_name%.ms}.clean.residual.tar.gz vos:MWSynthesis/VLA/archival/test_bad_cleans/channel_${channel}/

vcp ${TMPDIR}/proc_${channel}/${ms_name%.ms}.clean.flux.tar.gz vos:MWSynthesis/VLA/archival/test_bad_cleans/channel_${channel}/

vcp ${TMPDIR}/proc_${channel}/${ms_name%.ms}.clean.flux.pbcoverage.tar.gz vos:MWSynthesis/VLA/archival/test_bad_cleans/channel_${channel}/
