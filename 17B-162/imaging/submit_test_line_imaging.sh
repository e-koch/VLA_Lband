
# Pass variables to the job script to separately run each SPW

sbatch --account=rrg-eros-ab --export=spw='0' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --account=rrg-eros-ab --export=spw='1' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --account=rrg-eros-ab --export=spw='2' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --account=rrg-eros-ab --export=spw='3' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --account=rrg-eros-ab --export=spw='4' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --account=rrg-eros-ab --export=spw='5' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --account=rrg-eros-ab --export=spw='6' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --account=rrg-eros-ab --export=spw='7' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --account=rrg-eros-ab --export=spw='8' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --account=rrg-eros-ab --export=spw='9' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
