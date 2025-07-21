# APEX Monte-Carlo Generator
Generates Monte-Carlo data for the purpose of APEX optics calibration.

This project is a fork of both the Apex G4MC simulation package obtained from: (https://github.com/kvardan/APEX_G4MC), 
and also a fork of a modified version of the SIMC Hall-A simulation package: (https://github.com/JeffersonLab/simc_gfortran). 


### To Compile: 
Clone this repository into an empty driectory of your choice: 

    git clone https://github.com/sethcarl2000/apex_mc_generator

Go into both compiled directories, and use the macro to compile them. for G4MC:
  
    cd g4mc
    ./compile.sh -f
  
Here, the '-f' flag forces a rebuild of the build/ directory. Same for simc_hrs:
  
    cd simc_hrs
    ./compile.sh -f
  
Then, you must also uncompress the septum map file:
  
    gunzip g4mc/septum/Septa-JB_map.table.gz


## G4MC 
This is a fork of [this repo](https://github.com/kvardan/APEX_G4MC) by Vardan Khachatryan. G4MC simulates the propagation of electrons/positrons from the apex scattering chamber, through the septum, to the front of the Q1 quadrapole (the frist magnetic lens of the HRS). The simulation is built using the Geant4 simulation package. 

See the [script_gen-Q1-tracks](https://github.com/sethcarl2000/apex_mc_generator/edit/main/README.md#script_gen-q1-tracks) below to see how execution is performed. The following 'commands' must exist in the .mac file you provide when you invoke the G4MC executable. When you run the 'script_gen-Q1-tracks' script and specify a target to use (production or vertical wire), all of this is handled automatically for you. See below for usage of the 'script_gen-Q1-tracks' script. 

The name of the output .root file is specified via: 

```
/mydet/outfile output_file.root
```

Each time a new lepton is generated, its position, momentum and particle species (positron/electron) are defined in [HRSPrimaryGeneratorAction::GeneratePrimaries()](https://github.com/sethcarl2000/apex_mc_generator/blob/main/g4mc/src/HRSPrimaryGeneratorAction.cc). 
you can use the following commands to specify the 'zone' over which particles will be generated **Relative to the center of the APEX target scattering chamber, NOT the Hall-A centerline**. As an example, these vertex generation bounds are taken from [slurm_infiles/settings_target_production.txt](github.com/sethcarl2000/apex_mc_generator/blob/main/g4mc/slurm_infiles/settings_target_production.txt): 

```
/mydet/gunXLow  -3.5 mm
/mydet/gunXHigh  3.5 mm

/mydet/gunYLow  -3.5 mm
/mydet/gunYHigh  3.5 mm

/mydet/gunZLow  -240.0 mm
/mydet/gunZHigh  330.0 mm
```

And the following command will specify the range of lepton momenta to use: 

```
/mydet/particle1/momentumLow  1058 MeV
/mydet/particle1/momentumHigh 1149 MeV
```

To specify a range of particle trajectores, the particle picks a random spot on the face of the sieve to 'aim' at in the given range (for the right or left arm). The range is given relative to the central 'big hole' of either sieve plate: 

```
/mydet/sieveXLow   -85.0 mm
/mydet/sieveXHigh   85.0 mm

/mydet/sieveYLow   -80.0 mm
/mydet/sieveYHigh   70.0 mm
```

The HRS arm to use must also be specified: 

```
/mydet/use_RHRS true
```

Specifying the HRS arm to use will automatically generate positrons(electrons) and aim them at the Right(Left) HRS.  

The positions and angles of both sieve plates were surveyed during the APEX production run. This data is used to accurately set 'sieve-coordinates' relative to the Hall Centerline, and the APEX scattering chamber centerline. In [ApexTargetGeometry.hh](https://github.com/sethcarl2000/apex_mc_generator/blob/main/g4mc/src/ApexTargetGeometry.hh), you can see the 'ApexTargetGeometry' namespace, in which several hard-coded functions prodvide information about the position and rotation of the sieve relative to the APEX scattering chamber center, as well as the offset between the APEX target chamber and Hall Centerline. 

> [!NOTE]
> By default, each part of the G4MC simulation uses the **Hall Coordinate System** (HCS), in which all coordinates are given relative to the center of Hall-A. In this system, the z-axis points downbeam, the y-axis points upward normal to the floor, and the x-axis points left when you face downbeam. _However_, the target-chamber generation limits are defined relative to the APEX scattering chamber center. See the below section on the definition of output variables for more details on the different coordinate systems used.
> Unless otherwise specified, the G4MC simulation uses **mm for all displacements** and **MeV/c for all momenta** during the simulation process, but **output varaibles are in meters**. See below for more information. 



Rather than using sensitive detectors, this code was set up to check the 'physical volume' of a lepton at each 1mm simulation step. In words: for each 1mm step a lepton takes through the septum/scattering chamber, the code checks to see if the lepton has reached the physical volume defined for the Q1 front window. This is done in [HRSSteppingAction::UserSteppingAction()](https://github.com/sethcarl2000/apex_mc_generator/blob/main/g4mc/src/HRSSteppingAction.cc) 

### Output Variables

The output .root file will contain a TTree named **tracks_Q1**. The branches in this tree are described below. There will also be a TParameter\<bool\> object called **is_RHRS**. This is neeed by the simc_hrs code to determine which HRS arm configuration to simulate. 

Once a lepton reaches the Q1 front-window of either arm, the HRSSteppingAction::UserSteppingAction method defines the following set of output variables, via the [TFileHandler](https://github.com/sethcarl2000/apex_mc_generator/blob/main/g4mc/src/TFileHandler.cc) class: 

> [!NOTE]
> While the G4MC uses mm for displacements during the simulation process, all of these output variables use **meters for all displacements** and **MeV/c for all momenta** 

+ **position_vtx** (TVector3): This represents the position at which the particle was generated, <ins>Relative to the APEX scattering chamber center</ins>. Despite this shift from the Hall Centerline, the orientaiton of the _xyz_-axes are the same as in HCS: _z_-axis downbeam, _x_-axis pointing left as you face downbeam, and _y_-axis pointing up, normal to the floor.
+ **momentum_vtx** (TVector3): This represents the 3-momentum of the particle at its moment of generation. Axes are the HCS system described above.

+ **position_sieve** (TVector3): This vector is in the **Sieve Coordinate System** (SCS), which is a modified verison of the standard Target Coordinate System. SCS is defined in the following way:
    - The center of either sieve plane is the central 'big hole' (each sieve plate has 2 holes which are bigger than the rest). On the plane of the sieve which faces the scattering chamber, the origin of SCS is the center of the central big-hole's target-facing opening. Note that the position of this origin is hard-coded in the [ApexTargetGeometry namespace](https://github.com/sethcarl2000/apex_mc_generator/blob/e6167227437cf83579fa1241d5c4b3740539e1b4/g4mc/src/ApexTargetGeometry.hh), which uses data from the APEX sieve-survey. 
    - The _z_-axis is coliniear with the axis of the central big hole, facing from the target into the septum.
    - The _x_ and _y_ axes both lie in the plane of the sieve's face, however, **The _x_-axis points down, normal to the hall floor, and the _y_-axis points to the left, when facing the sieve from the target.**
+ **momentum_sieve** (TVector3): Oriented using the SCS axes described above.

+ **position_Q1** (TVector3): This vector is in the "traditional" **Target Coordinate System** (TCS).
    - The origin is the Hall Centerline.
    - As in the Sieve Coordinate System, The _x_-axis points normal-down to the hall floor, and the _y_-axis points to the left, when facing the spectrometer opening from the target. Note that while each sieve plane _z_-axis is rotated from the nominal HCS _z_-axis at an agle of ~ +/- 5 deg, the HRS arms in the APEX configuration are at angles of +/- 12.5 deg. Therefore, the _z_- and _y_-axes of the SCS and TCS are distinct. 
+ **momentum_Q1** (TVector3): Oriented using the TCS described above. 


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



  

  
