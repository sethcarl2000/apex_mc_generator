#!/bin/bash

###########################################
#
#  Launches a test-run of the SIMC simulation of the HRS.
#
#
#  To convert a Q1_front.dat txt file to SIMC-ready output, invoke: 
#
#  > root -l -b -q 'convert_Q1_to_simc.C'
#
#  then call this script. 
#
###########################################

if [ "$1" = "-r" ]
then
    echo -n "Converting 'Q1_front.dat' txt-file to simc-ready data..."

    root -l -b -q 'convert_Q1_to_simc.C' >/dev/null

    echo "done."
    echo
fi

echo "Running SIMC..." 

cd $WORK_MC/simc
./mc_hrs_single < inp.txt
