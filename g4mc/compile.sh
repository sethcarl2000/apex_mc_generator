#!/bin/bash

cd $WORK_MC/g4mc

if [ "$1" = "-f" ]
then
    echo "deleting & remaking directory 'build'..."
    rm -rf build
    mkdir build
fi
    
#rm -rf build
#mkdir build

cmake -B build -S . 
cmake --build build -j8 
