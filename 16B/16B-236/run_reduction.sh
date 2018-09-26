
# Reduce all of the line SPWs from 16B-236
# Running 4 tracks at a time


two36_path=$HOME/space/ekoch/VLA_tracks/16B-236/reduction

cd $two36_path

pids=

echo 16B-236_11_11_16
cd 16B-236_11_11_16
(tar -xf 16B-236.sb32623829.eb32998154.57703.975014375.tar && 16B-236.sb32623829.eb32998154.57703.975014375.tar && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-236.sb32623829.eb32998154.57703.975014375 F lines && cd $two36_path/16B-236_11_11_16/16B-236_11_11_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-236.sb32623829.eb32998154.57703.975014375.speclines.ms) &
pids+=" $!"


echo 16B-236_11_16_16
cd $two36_path/16B-236_11_16_16
(tar -xf 16B-236.sb32623829.eb33005682.57708.9880922338.tar && rm 16B-236.sb32623829.eb33005682.57708.9880922338.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-236.sb32623829.eb33005682.57708.9880922338 F lines && cd $two36_path/16B-236_11_16_16/16B-236_11_16_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-236.sb32623829.eb33005682.57708.9880922338.speclines.ms) &
pids+=" $!"

echo 16B-236_11_26_16
cd $two36_path/16B-236_11_26_16
(tar -xf 16B-236.sb32623829.eb33039582.57718.93170724537.tar && rm 16B-236.sb32623829.eb33039582.57718.93170724537.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-236.sb32623829.eb33039582.57718.93170724537 F lines && cd $two36_path/16B-236_11_26_16/16B-236_11_26_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-236.sb32623829.eb33039582.57718.93170724537.speclines.ms) &
pids+=" $!"

echo 16B-236_11_27_16
cd $two36_path/16B-236_11_27_16
(tar -xf 16B-236.sb32623829.eb33039653.57719.98023209491.tar && rm 16B-236.sb32623829.eb33039653.57719.98023209491.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-236.sb32623829.eb33039653.57719.98023209491 F lines && cd $two36_path/16B-236_11_27_16/16B-236_11_27_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-236.sb32623829.eb33039653.57719.98023209491.speclines.ms) &
pids+=" $!"

wait $pids || { echo "There was an error" >&2; exit 1; }

echo "All CASA jobs exited. 1/3"

pids=

echo 16B-236_11_28_16
cd $two36_path/16B-236_11_28_16
(tar -xf 16B-236.sb32623829.eb33043297.57720.92533724537.tar && rm 16B-236.sb32623829.eb33043297.57720.92533724537.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-236.sb32623829.eb33043297.57720.92533724537 F lines && cd $two36_path/16B-236_11_28_16/16B-236_11_28_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-236.sb32623829.eb33043297.57720.92533724537.speclines.ms) &
pids+=" $!"

echo 16B-236_11_29_16
cd $two36_path/16B-236_11_29_16
(tar -xf 16B-236.sb32658813.eb33043300.57721.2355058449.tar && rm 16B-236.sb32658813.eb33043300.57721.2355058449.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-236.sb32658813.eb33043300.57721.2355058449 F lines && cd $two36_path/16B-236_11_29_16/16B-236_11_29_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-236.sb32658813.eb33043300.57721.2355058449.speclines.ms) &
pids+=" $!"

echo 16B-236_12_08_16
cd $two36_path/16B-236_12_08_16
(tar -xf 16B-236.sb32623829.eb33067127.57730.91522394676.tar && rm 16B-236.sb32623829.eb33067127.57730.91522394676.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-236.sb32623829.eb33067127.57730.91522394676 F lines && cd $two36_path/16B-236_12_08_16/16B-236_12_08_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-236.sb32623829.eb33067127.57730.91522394676.speclines.ms) &
pids+=" $!"

echo 16B-236_12_09_16
cd $two36_path/16B-236_12_09_16
(tar -xf 16B-236.sb32623829.eb33069261.57731.904499571756.tar && rm 16B-236.sb32623829.eb33069261.57731.904499571756.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-236.sb32623829.eb33069261.57731.904499571756 F lines && cd $two36_path/16B-236_12_09_16/16B-236_12_09_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-236.sb32623829.eb33069261.57731.904499571756.speclines.ms) &
pids+=" $!"

wait $pids || { echo "There was an error" >&2; exit 1; }

echo "All CASA jobs exited. 2/3"

pids=

echo 16B-236_12_15_16
cd $two36_path/16B-236_12_15_16
(tar -xf 16B-236.sb32623829.eb33077539.57737.926878715276.tar && rm 16B-236.sb32623829.eb33077539.57737.926878715276.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-236.sb32623829.eb33077539.57737.926878715276 F lines && cd $two36_path/16B-236_12_15_16/16B-236_12_15_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-236.sb32623829.eb33077539.57737.926878715276.speclines.ms) &
pids+=" $!"

echo 16B-236_12_17_16
cd $two36_path/16B-236_12_17_16
(tar -xf 16B-236.sb32658813.eb33081823.57739.161825104165.tar && rm 16B-236.sb32658813.eb33081823.57739.161825104165.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-236.sb32658813.eb33081823.57739.161825104165 F lines && cd $two36_path/16B-236_12_17_16/16B-236_12_17_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-236.sb32658813.eb33081823.57739.161825104165.speclines.ms) &
pids+=" $!"

echo 16B-236_12_30_16
cd $two36_path/16B-236_12_30_16
(tar -xf 16B-236.sb33078950.eb33143443.57752.87391505787.tar && rm 16B-236.sb33078950.eb33143443.57752.87391505787.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-236.sb33078950.eb33143443.57752.87391505787 F lines && cd $two36_path/16B-236_12_30_16/16B-236_12_30_16_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-236.sb33078950.eb33143443.57752.87391505787.speclines.ms) &
pids+=" $!"

echo 16B-236_01_01_17
cd $two36_path/16B-236_01_01_17
(tar -xf 16B-236.sb32658813.eb33143793.57754.172903252314.tar && rm 16B-236.sb32658813.eb33143793.57754.172903252314.tar  && ~/casa-release-5.1.2-4.el7/bin/casa -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/ms_split.py 16B-236.sb32658813.eb33143793.57754.172903252314 F lines && cd $two36_path/16B-236_01_01_17/16B-236_01_01_17_speclines &&  ~/casa-release-5.1.2-4.el7/bin/casa --pipeline -c ~/Dropbox/code_development/VLA_Lband/16B/pipeline_scripts/casa_pipeline_lines.py 16B-236.sb32658813.eb33143793.57754.172903252314.speclines.ms) &
pids+=" $!"

wait $pids || { echo "There was an error" >&2; exit 1; }

echo "All CASA jobs exited. 3/3"
