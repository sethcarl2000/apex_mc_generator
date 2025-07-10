# APEX_monte_carlo
Generates Monte-Carlo data for the purpose of APEX optics calibration.

**NOTE: CODE WILL RUN IMPROPERLY IF NOT DONE**

The following file was too large to upload to github, so I had to compress it: 

> g4mc/septum/Septa-JB_map.table 

if you clone this repo and want to execute it, you must first run the following command:

    gunzip g4mc/septum/Septa-JB_map.table.gz

You don't have to recompile anything, it should work fine. 


This project is a fork of both the Apex G4MC simulation package obtained from: (https://github.com/kvardan/APEX_G4MC), 
and also a fork of a modified version of the SIMC Hall-A simulation package: (https://github.com/JeffersonLab/simc_gfortran). 

Usage / code: 


## script_gen-Q1-tracks 

This script is a self-contained example of how to execute the G4MC simulation of the APEX scattering chamber and Septum magnet. 

### Inputs/Arguments: 
There are three possible positional arguments for this script, the first two are necessary; the script will refuse to execute if either of the first
two are missing. 

- [R/L] arm to use; either RHRS or LHRS. 

- [_number of monte-carlo events_] a positive integer which is the number of events to simulate. 

- [_target_] (optional) The target to use; as of writing, valid inputs are 'proudction','V1','V2','V3'. 'production' is default if none is specified.   
  This parameter changes the range over which the paritcles are generated. For 'production', particles are uniformly generated in a 'box' which correpsonds 
  to the approximate extent of the APEX production target, which was 10 tungsten foils spaced evenly along the z-axis of the scattering chamber. For each 
  of the vertical wire targets (V1,V2,V3), the generation of particles is restricted to the survey-verified positions of the vertical wires. particles are
  still generated in a 2mm uniform range in y_HCS to account for the raster used for optics-wire runs. 

See example below to see examples of proper usage. 

### Ouputs:
There are two paths hard-coded in; feel free to change these to suit your needs.* 

    VOL_DIR="/volatile/halla/apex/mc" 
    PATH_G4MC="/work/halla/apex/disk1/sethhall/apex_mc_generator/g4mc"

**Note**: PATH_G4MC must point to the g4mc source directory, this script will fail execution otherwise! 

Upon successful execution of the script, a new directory will be created:  

    ${VOL_DIR}/temp/job_${run_id}
    
**Note**: if this directory already exists, all of its contents will be overwritten.     

The output root file will be placed under: 

    ${VOL_DIR}/tracks_Q1/out_${hrs_arm}_${target}_${run_id}.root

Where hrs_arm="RHRS"/"LHRS", and run_id is '0000' if the script is executed on the command line; otherwise, it is the slurm JobID. 
If the script is submitted to run as a slurm job, the out/err logs will be under the directory: 

    ${VOL_DIR}/temp/slurm

### Execution Examples: 
All of these are executed in a bash shell on ifarm: 

- This simulates 100k positrons, aimed at the RHRS-side of the septum magnet (the target used is production): 

      ./script_gen-Q1-tracks R 100000

  This will generate a new directory: 

      ${VOL_DIR}/temp/job_0000
    
  And an associated output root file:

      ${VOL_DIR}/tracks_Q1/out_RHRS_R_production_0000.root

- This simulates 250k electrons, generating them at the position of the V1 vertical optics wire: 

      ./scirpt_gen-Q1-tracks L 250000 V1

  This will generate a new directory: 

      ${VOL_DIR}/temp/job_0000
    
  And an associated output root file:

      ${VOL_DIR}/tracks_Q1/out_LHRS_V1_0000.root
    

- This _tests_ the submission of the script for 200k events (LHRS), with 90 mins of runtime. the '--test-only --partition=production' parameter
  means that this job won't actually submit; instead, it will tell you when the job would _start_ running, if you submitted it at this moment, with these 
  parameters, on the 'production' partition: 

      sbatch --test-only --partition=production --job-name=apex_septum_sim_L --time=90 script_gen-Q1-tracks L 200000

- This would actually submit a job to execute this script: 

      sbatch --job-name=apex_septum_sim_R_V2 --time=100 script_gen-Q1-tracks R 200000 V2

  The slurm job Id will be an 8-digit number like 50501919, and thus the new temp-file directory will be: 

      ${VOL_DIR}/temp/job_50501919
    
  And associated ouput root file:

      ${VOL_DIR}/tracks_Q1/out_RHRS_V2_50501919.root

**NOTE**: you must specify the '--time' option, as no time option is specified by default in the script. To estimate the amount of time needed to run, I use
the fact that the script takes about 17 ms/event for its simulation; I usually round up to 25 ms/event to leave some buffer room. _If the job runs over the
time budget you give it, it will fail!_ 



  

  