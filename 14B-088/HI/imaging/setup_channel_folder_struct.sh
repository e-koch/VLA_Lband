
# Create separate directory structure for each channel
start=670
nchan=400

cd /home/ekoch/m33/14B-088/single_channels/

for i in {start..$((start+nchan))}; do
    mkdir channel_${i}
    mv ../channel_ms/14B-088_HI_LSRK.ms.contsub_channel_${i}.ms .
    mv ../mask_channels/M33_14B-088_HI_mask_channel_${i}.image .
    mv ../model_channels/M33_14B-088_HI_model_channel_${i}.image .
