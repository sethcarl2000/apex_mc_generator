#!/bin/bash


array_size=10
n_events_per_array=50000000

for ARM in LHRS RHRS;
do
    for WIRE in V1 V2 V3;
    do

	cmd="sbatch --job-name=simulate_sieve_${ARM}_${WIRE} --time=600 --array=1-$array_size slurm_script ${ARM} ${WIRE} all_types $n_events_per_array"

	echo "cmd: ${cmd}"

	eval "$cmd"
			
    done
done
