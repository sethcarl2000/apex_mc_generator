#!/bin/bash

if [ "$1" = "-f" ]
then
    rm -rf build
    mkdir build
fi

cmake -B build -S .
cmake --build build

