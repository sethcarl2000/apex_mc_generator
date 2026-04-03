#!/bin/bash

STARTING_THICKNESS=${1}

THICKNESS_SPACING=${2}

N_RUNS=${3}

N_EVENTS=200000

echo "Executing ${3} runs, starting with a thickness of ${1} um, and each run spaced by ${2} um"

SIZE_UM=${STARTING_THICKNESS}

echo "starting size: ${SIZE_MM} mm"

i=0
while [[ $i < ${N_RUNS} ]]; do

    echo "size: ${SIZE_MM} mm" 

    PATH_MACRO="myrun_${SIZE_UM}um.mac"
    cp run_variable_thickness.mac ${PATH_MACRO}

    PATH_OUTFILE="data/test_target_${SIZE_UM}um.root"

    echo "" >> ${PATH_MACRO}
    echo "# output file " >> ${PATH_MACRO}
    echo "/global/path_outfile ${PATH_OUTFILE}" >> ${PATH_MACRO}

    SIZE_MM=$(echo "scale=3; ${SIZE_UM}/1000" | bc)
    echo "" >> ${PATH_MACRO}
    echo "# target thickness" >> ${PATH_MACRO}
    echo "/detector/target_thickness ${SIZE_MM} mm" >> ${PATH_MACRO}
    
    echo "" >> ${PATH_MACRO}
    echo "# Start the run" >> ${PATH_MACRO}
    echo "/run/initialize" >> ${PATH_MACRO}
    echo "/run/beamOn ${N_EVENTS}" >> ${PATH_MACRO}

    time ./exampleB1 ${PATH_MACRO} > out_${SIZE_UM}.log

    hadd -j 4 "${PATH_OUTFILE}" data/*um_t*
    rm data/*um_t*

    i=$(( i + 1 ))        
    SIZE_UM=$(( SIZE_UM + THICKNESS_SPACING ))
    
done
