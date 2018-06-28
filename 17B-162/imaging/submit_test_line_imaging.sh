
# Pass variables to the job script to separately run each SPW

sbatch --export=spw='0' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --export=spw='1' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --export=spw='2' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --export=spw='3' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --export=spw='4' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --export=spw='5' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --export=spw='6' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --export=spw='7' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --export=spw='8' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
sbatch --export=spw='9' $HOME/code/VLA_Lband/17B-162/imaging/test_line_imaging.sh
