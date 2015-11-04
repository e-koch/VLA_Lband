
## OH Lines

cd /media/eric/Data_1/14B-088_OH/OH1612/

for ms in *.ms; do
    casa-4.2.2 -c /home/eric/Dropbox/code_development/VLA_Lband/14B-088/Lines/plot_scans.py $ms M33* 0
done

mkdir scan_plots
mv *_scan_plots scan_plots
rm *.log

cd /media/eric/Data_1/14B-088_OH/OH1665/

for ms in *.ms; do
    casa-4.2.2 -c /home/eric/Dropbox/code_development/VLA_Lband/14B-088/Lines/plot_scans.py $ms M33* 0
done

mkdir scan_plots
mv *_scan_plots scan_plots
rm *.log

cd /media/eric/Data_1/14B-088_OH/OH1667/

for ms in *.ms; do
    casa-4.2.2 -c /home/eric/Dropbox/code_development/VLA_Lband/14B-088/Lines/plot_scans.py $ms M33* 0
done

mkdir scan_plots
mv *_scan_plots scan_plots
rm *.log

cd /media/eric/Data_1/14B-088_OH/OH1720/

for ms in *.ms; do
    casa-4.2.2 -c /home/eric/Dropbox/code_development/VLA_Lband/14B-088/Lines/plot_scans.py $ms M33* 0
done

mkdir scan_plots
mv *_scan_plots scan_plots
rm *.log

## HRLs

cd /media/eric/Data_1/14B-088_HRL/H152alp/

for ms in *.ms; do
    casa-4.2.2 -c /home/eric/Dropbox/code_development/VLA_Lband/14B-088/Lines/plot_scans.py $ms M33* 0
done

mkdir scan_plots
mv *_scan_plots scan_plots
rm *.log

cd /media/eric/Data_1/14B-088_HRL/H153alp/

for ms in *.ms; do
    casa-4.2.2 -c /home/eric/Dropbox/code_development/VLA_Lband/14B-088/Lines/plot_scans.py $ms M33* 0
done

mkdir scan_plots
mv *_scan_plots scan_plots
rm *.log

cd /media/eric/Data_1/14B-088_HRL/H158alp/

for ms in *.ms; do
    casa-4.2.2 -c /home/eric/Dropbox/code_development/VLA_Lband/14B-088/Lines/plot_scans.py $ms M33* 0
done

mkdir scan_plots
mv *_scan_plots scan_plots
rm *.log

cd /media/eric/Data_1/14B-088_HRL/H164alp/

for ms in *.ms; do
    casa-4.2.2 -c /home/eric/Dropbox/code_development/VLA_Lband/14B-088/Lines/plot_scans.py $ms M33* 0
done

mkdir scan_plots
mv *_scan_plots scan_plots
rm *.log

cd /media/eric/Data_1/14B-088_HRL/H166alp/

for ms in *.ms; do
    casa-4.2.2 -c /home/eric/Dropbox/code_development/VLA_Lband/14B-088/Lines/plot_scans.py $ms M33* 0
done

mkdir scan_plots
mv *_scan_plots scan_plots
rm *.log
