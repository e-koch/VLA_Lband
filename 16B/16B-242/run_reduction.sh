
# Reduce all of the line SPWs from 16B-242
# Running 4 tracks at a time


two42_path=$HOME/space/ekoch/VLA_tracks/16B-242/reduction
two42_path=$HOME/bigdata/ekoch/VLA_tracks/16B-242/reduction

cd $two42_path

pids=

echo 16B-242_10_12_16
cd 16B-242_10_12_16
(tar -xf 16B-242.sb32681213.eb32949259.57673.06817918981.tar && rm 16B-242.sb32681213.eb32949259.57673.06817918981.tar && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-242.sb32681213.eb32949259.57673.06817918981 F lines && cd $two42_path/16B-242_10_12_16/16B-242_10_12_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-242.sb32681213.eb32949259.57673.06817918981.speclines.ms) &
pids+=" $!"


echo 16B-242_10_13_16
cd 16B-242_10_13_16
(tar -xf 16B-242.sb32681213.eb32954876.57674.07753513889.tar && rm 16B-242.sb32681213.eb32954876.57674.07753513889.tar && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-242.sb32681213.eb32954876.57674.07753513889 F lines && cd $two42_path/16B-242_10_13_16/16B-242_10_13_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-242.sb32681213.eb32954876.57674.07753513889.speclines.ms) &
pids+=" $!"

echo 16B-242_10_14_16
cd 16B-242_10_14_16
(tar -xf 16B-242.sb32681213.eb32956712.57675.0472325463.tar && rm 16B-242.sb32681213.eb32956712.57675.0472325463.tar && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-242.sb32681213.eb32956712.57675.0472325463 F lines && cd $two42_path/16B-242_10_14_16/16B-242_10_14_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-242.sb32681213.eb32956712.57675.0472325463.speclines.ms) &
pids+=" $!"

echo 16B-242_10_16_16
cd $two42_path/16B-242_10_16_16
(tar -xf 16B-242.sb32681213.eb32958324.57677.085873668984.tar && rm 16B-242.sb32681213.eb32958324.57677.085873668984.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-242.sb32681213.eb32958324.57677.085873668984 F lines && cd $two42_path/16B-242_10_16_16/16B-242_10_16_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-242.sb32681213.eb32958324.57677.085873668984.speclines.ms) &
pids+=" $!"

wait $pids || { echo "There was an error" >&2; exit 1; }

echo "All CASA jobs exited. 1/3"

pids=

echo 16B-242_10_20_16
cd $two42_path/16B-242_10_20_16
(tar -xf 16B-242.sb32681213.eb32961984.57681.03943858796.tar && rm 16B-242.sb32681213.eb32961984.57681.03943858796.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-242.sb32681213.eb32961984.57681.03943858796 F lines && cd $two42_path/16B-242_10_20_16/16B-242_10_20_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-242.sb32681213.eb32961984.57681.03943858796.speclines.ms) &
pids+=" $!"

echo 16B-242_10_25_16
cd $two42_path/16B-242_10_25_16
(tar -xf 16B-242.sb32681213.eb32969161.57686.01854777778.tar && rm 16B-242.sb32681213.eb32969161.57686.01854777778.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-242.sb32681213.eb32969161.57686.01854777778 F lines && cd $two42_path/16B-242_10_25_16/16B-242_10_25_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-242.sb32681213.eb32969161.57686.01854777778.speclines.ms) &
pids+=" $!"

# The remaining tracks were taken after M33 transit zenith. The obs setup used the 3C48 scan 2 for
# phase and amp cal. Optionally adjust the scan intents are untar-ing the file to avoid this.

echo 16B-242_11_01_16
cd $two42_path/16B-242_11_01_16
(tar -xf 16B-242.sb32614458.eb32978118.57693.289744062495.tar && rm 16B-242.sb32614458.eb32978118.57693.289744062495.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-242.sb32614458.eb32978118.57693.289744062495 F lines && cd $two42_path/16B-242_11_01_16/16B-242_11_01_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-242.sb32614458.eb32978118.57693.289744062495.speclines.ms) &
pids+=" $!"

echo 16B-242_11_02_16
cd $two42_path/16B-242_11_02_16
(tar -xf 16B-242.sb32614458.eb32978942.57694.29761952546.tar && rm 16B-242.sb32614458.eb32978942.57694.29761952546.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-242.sb32614458.eb32978942.57694.29761952546 F lines && cd $two42_path/16B-242_11_02_16/16B-242_11_02_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-242.sb32614458.eb32978942.57694.29761952546.speclines.ms) &
pids+=" $!"

wait $pids || { echo "There was an error" >&2; exit 1; }

echo "All CASA jobs exited. 2/3"

pids=

echo 16B-242_11_03_16
cd $two42_path/16B-242_11_03_16
(tar -xf 16B-242.sb32614458.eb32980749.57695.35484743056.tar && rm 16B-242.sb32614458.eb32980749.57695.35484743056.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-242.sb32614458.eb32980749.57695.35484743056 F lines && cd $two42_path/16B-242_11_03_16/16B-242_11_03_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-242.sb32614458.eb32980749.57695.35484743056.speclines.ms) &
pids+=" $!"

echo 16B-242_11_04_16
cd $two42_path/16B-242_11_04_16
(tar -xf 16B-242.sb32614458.eb32982501.57696.352311215276.tar && rm 16B-242.sb32614458.eb32982501.57696.352311215276.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-242.sb32614458.eb32982501.57696.352311215276 F lines && cd $two42_path/16B-242_11_04_16/16B-242_11_04_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-242.sb32614458.eb32982501.57696.352311215276.speclines.ms) &
pids+=" $!"

echo 16B-242_11_05_16
cd $two42_path/16B-242_11_05_16
(tar -xf 16B-242.sb32614458.eb32984320.57697.291263148145.tar && rm 16B-242.sb32614458.eb32984320.57697.291263148145.tar && tar -xf 16B-242_sb32614458_5.57697.30203707176.tar && rm 16B-242_sb32614458_5.57697.30203707176.tar && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/ms_split_16B-242_11_05.py F lines && cd $two42_path/16B-242_11_05_16/16B-242_11_05_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-242.sb32614458.eb32984320.57697.291263148145.speclines.ms) &
pids+=" $!"

echo 16B-242_11_06_16
cd $two42_path/16B-242_11_06_16
(tar -xf 16B-242.sb32614458.eb32984330.57698.29168733796.tar && rm 16B-242.sb32614458.eb32984330.57698.29168733796.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-242.sb32614458.eb32984330.57698.29168733796 F lines && cd $two42_path/16B-242_11_06_16/16B-242_11_06_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-242.sb32614458.eb32984330.57698.29168733796.speclines.ms) &
pids+=" $!"

wait $pids || { echo "There was an error" >&2; exit 1; }

echo "All CASA jobs exited. 3/3"

# Continuum reduction

export repo_path='/home/eric/ownCloud/code_development/VLA_Lband/'
export repo_path='/home/ekoch/ownCloud/code_development/VLA_Lband/'


export track_folder='16B-242_10_12_16'
export track_name='16B-242.sb32681213.eb32949259.57673.06817918981'

# Move the SDM file into the track folder
mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/16B/pipeline_scripts/ms_split.py $track_name F cont
# Clean up intermediate data
# rm -f "${track_name}.tar"
# rm -rf "${track_name}"
cd "${track_folder}_continuum"
# And reduce with the 5.4.1 pipeline
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/16B/pipeline_scripts/casa_pipeline_continuum_v541.py "${track_name}.continuum.ms"
cd ../../

export track_folder='16B-242_10_13_16'
export track_name='16B-242.sb32681213.eb32954876.57674.07753513889'

# Move the SDM file into the track folder
mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/16B/pipeline_scripts/ms_split.py $track_name F cont
# Clean up intermediate data
# rm -f "${track_name}.tar"
# rm -rf "${track_name}"
cd "${track_folder}_continuum"
# And reduce with the 5.4.1 pipeline
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/16B/pipeline_scripts/casa_pipeline_continuum_v541.py "${track_name}.continuum.ms"
cd ../../

export track_folder='16B-242_10_14_16'
export track_name='16B-242.sb32681213.eb32956712.57675.0472325463'

# Move the SDM file into the track folder
mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/16B/pipeline_scripts/ms_split.py $track_name F cont
# Clean up intermediate data
# rm -f "${track_name}.tar"
# rm -rf "${track_name}"
cd "${track_folder}_continuum"
# And reduce with the 5.4.1 pipeline
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/16B/pipeline_scripts/casa_pipeline_continuum_v541.py "${track_name}.continuum.ms"
cd ../../

export track_folder='16B-242_10_16_16'
export track_name='16B-242.sb32681213.eb32958324.57677.085873668984'

# Move the SDM file into the track folder
mv "${track_name}.tar" $track_folder
cd $track_folder
tar -xf "${track_name}".tar
~/casa-release-5.4.1-32.el7/bin/casa --nogui --log2term -c $repo_path/16B/pipeline_scripts/ms_split.py $track_name F cont
# Clean up intermediate data
# rm -f "${track_name}.tar"
# rm -rf "${track_name}"
cd "${track_folder}_continuum"
# And reduce with the 5.4.1 pipeline
~/casa-release-5.4.1-32.el7/bin/casa --pipeline --nogui --log2term -c $repo_path/16B/pipeline_scripts/casa_pipeline_continuum_v541.py "${track_name}.continuum.ms"
cd ../../
