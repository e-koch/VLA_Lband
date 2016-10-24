
# Create the symbolic links for my custom pipeline version to work.
# This gives access to the general imaging scripts and imaging utility functions.

# You should also add the new path to paths.py

echo "Run this script in the new pipeline directory"

# ln -s $HOME/Dropbox/code_development/VLA_Lband/paths.py
# ln -s $HOME/Dropbox/code_development/VLA_Lband/CASA_functions
# ln -s $HOME/Dropbox/code_development/VLA_Lband/imaging_pipeline

ln -s ../../paths.py
ln -s ../../CASA_functions
ln -s ../../imaging_pipeline
ln -s ../../flagging_scripts
