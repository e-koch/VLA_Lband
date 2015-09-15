#!/bin/bash

# Cleans-up a given directory with the outputs of the EVLA pipeline (v1.3.0)

pipe_direc=$1
data_type=$2

cd $pipe_direc

rm -rf *.last

tar -zcf calibrators.ms.tar.gz calibrators.ms

mv *.log logs/

tar -zcf logs.tar.gz logs

rm -rf BPcal.b BPinitialgain.g finalBPinitialgain.g finaldelayinitialgain.g fluxphaseshortgaincal.g phaseshortgaincal.g
rm -rf semiFinaldelayinitialgain.g testBPcal.b testBPdinitialgain.g testdelayinitialgain.g testdelay.k testgaincal.g

mv scan_plots weblog

tar -zcf weblog.$data_type.tar.gz weblog

tar -zcf test_images.tar.gz test_images