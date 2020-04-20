
# Run pipeline on continuum SPWs

export repo_path='/home/ekoch/ownCloud/code_development/VLA_Lband/'

export track_folder='17B-162_09_23_17'
export track_name='17B-162.sb34051874.eb34497862.58019.12101613426'

# Move the SDM file into the track folder
mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F cont
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_continuum"
# And reduce with the 5.4.1 pipeline
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/casa_pipeline_continuum.py "${track_name}.continuum.ms"
cd ../../

export track_folder='17B-162_09_26_17'
export track_name='17B-162.sb34051874.eb34512512.58022.09389369213'

# Move the SDM file into the track folder
mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F cont
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_continuum"
# And reduce with the 5.4.1 pipeline
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/casa_pipeline_continuum.py "${track_name}.continuum.ms"
cd ../../