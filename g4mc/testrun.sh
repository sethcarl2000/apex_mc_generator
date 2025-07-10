#!/bin/bash

cd /work/halla/apex/disk1/sethhall/apex_mc_generator/g4mc

n_trials=100


# reset all contents of /run/
rm -rf run/*

cp -a files_in/* run/.

ln -s $WORK_MC/g4mc/build/G4MC run/G4MC

cd run

#see if the user inputted a custom number of trials
if [ "$1" != "" ]
then
    rm beam.mac
    
    n_trials=$1

    echo "Executing $n_trials...." 
    
    echo "/run/beamOn $n_trials" > beam.mac
fi    


./G4MC -m 1 run_settings.mac beam.mac

ls -lrth

