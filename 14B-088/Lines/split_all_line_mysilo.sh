
# Assume I have my data storage on Silo mounted


# RRLs. Skip the one at ~1.2 GHz. It sits on a lovely time-variable RFI peak...
echo "RRLs"
echo "1/9"
casa -c ~/Dropbox/code_development/VLA_Lband/14B-088/Lines/split_spw_from_all_tracks.py 3 target /home/eric/Jasper/ /media/eric/Data_1/14B-088_HRL/
echo "2/9"
casa -c ~/Dropbox/code_development/VLA_Lband/14B-088/Lines/split_spw_from_all_tracks.py 4 target /home/eric/Jasper/ /media/eric/Data_1/14B-088_HRL/
echo "3/9"
casa -c ~/Dropbox/code_development/VLA_Lband/14B-088/Lines/split_spw_from_all_tracks.py 6 target /home/eric/Jasper/ /media/eric/Data_1/14B-088_HRL/
echo "4/9"
casa -c ~/Dropbox/code_development/VLA_Lband/14B-088/Lines/split_spw_from_all_tracks.py 10 target /home/eric/Jasper/ /media/eric/Data_1/14B-088_HRL/
echo "5/9"
casa -c ~/Dropbox/code_development/VLA_Lband/14B-088/Lines/split_spw_from_all_tracks.py 11 target /home/eric/Jasper/ /media/eric/Data_1/14B-088_HRL/

# OH
echo "OH"
echo "6/9"
casa -c ~/Dropbox/code_development/VLA_Lband/14B-088/Lines/split_spw_from_all_tracks.py 5 target /home/eric/Jasper/ /media/eric/Data_1/14B-088_OH/
echo "7/9"
casa -c ~/Dropbox/code_development/VLA_Lband/14B-088/Lines/split_spw_from_all_tracks.py 7 target /home/eric/Jasper/ /media/eric/Data_1/14B-088_OH/
echo "8/9"
casa -c ~/Dropbox/code_development/VLA_Lband/14B-088/Lines/split_spw_from_all_tracks.py 8 target /home/eric/Jasper/ /media/eric/Data_1/14B-088_OH/
echo "9/9"
casa -c ~/Dropbox/code_development/VLA_Lband/14B-088/Lines/split_spw_from_all_tracks.py 9 target /home/eric/Jasper/ /media/eric/Data_1/14B-088_OH/
