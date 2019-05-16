
# Restore the data given the pipeline products
# Expect the folder structure to already be there

export repo_path='/home/ekoch/ownCloud/code_development/VLA_Lband/'
export backup_path='/home/ekoch/work2/ekoch/VLA_backups/17B-162/Lines/'


export track_folder='17B-162_09_23_17'
export track_name='17B-162.sb34051874.eb34497862.58019.12101613426'

# Move the SDM file into the track folder
mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_09_26_17'
export track_name='17B-162.sb34051874.eb34512512.58022.09389369213'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_10_16_17'
export track_name='17B-162.sb34293636.eb34601252.58042.33489737268'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_10_17_17_early'
export track_name='17B-162.sb34051874.eb34603795.58043.05213521991'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_10_17_17_late'
export track_name='17B-162.sb34293636.eb34603801.58043.328202476856'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_10_19_17'
export track_name='17B-162.sb34293636.eb34618017.58045.33243515046'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_10_21_17'
export track_name='17B-162.sb34293636.eb34620280.58047.35365534722'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_10_22_17_early'
export track_name='17B-162.sb34051874.eb34620290.58048.0411411574'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_10_22_17_late'
export track_name='17B-162.sb34293636.eb34620658.58048.32964618056'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_10_26_17'
export track_name='17B-162.sb34293636.eb34630839.58052.310646990736'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_10_28_17'
export track_name='17B-162.sb34051874.eb34635762.58054.009438877314'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_10_29_17'
export track_name='17B-162.sb34293636.eb34635780.58055.299382048615'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_10_30_17'
export track_name='17B-162.sb34293636.eb34635792.58056.29234890046'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_10_31_17'
export track_name='17B-162.sb34293636.eb34642622.58057.30859833333'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_11_19_17'
export track_name='17B-162.sb34051874.eb34701346.58076.96044393518'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_11_25_17'
export track_name='17B-162.sb34051874.eb34742796.58082.91476392361'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../

export track_folder='17B-162_12_01_17'
export track_name='17B-162.sb34051874.eb34756401.58088.93220465278'

mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/17B-162/pipeline_scripts/ms_split.py $track_name F lines
# Clean up intermediate data
rm -f "${track_name}.tar"
rm -rf "${track_name}.ms"*
cd "${track_folder}_speclines"
# Copy over pipeline products.
# This grabs the flags and the caltables
cp "${backup_path}/${track_folder}"/*.tgz .
# This grabs the applycal command to use
cp "${backup_path}/${track_folder}"/*.txt .
# And restore with the 5.4.1 pipeline!
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/vla_pipeline_restore.py "${track_name}.speclines.ms"
cd ../../
