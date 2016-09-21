
# Run imaging tests w/ different parameters/CASA versions

filename=~/Dropbox/code_development/VLA_Lband/14B-088/HI/imaging/imaging_tests/HI_testing_channel_clean.py

ms_name=
modelname=
maskname=

declare -a versions=("casa-4.2.2" "casa-4.3.1" "casa-4.4" "casa-4.5" "casa")

# Number of runs to do at once
N=4


## now loop through the above array
for version in "${versions[@]}"; do

    echo $version

    for model in None modelname; do

        for mask in None maskname; do

            for use_fields in T F; do

                for use_mscale in T F; do

                    for use_tclean in T F; do

                        $version -c $filename ms_name model mask $use_fields $use_mscale $use_tclean

                    done

                done

            done

        done

    done

done
# Start combining all combination of parameters
casa-4.2.2 --logfile XXX -c $filename XXX