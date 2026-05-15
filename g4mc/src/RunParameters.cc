
#include "RunParameters.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exception.hh"
#include "TString.h"

#include <map> 

//______________________________________________________________________________
RunParameters::RunParameters()
{
    //add the commands we're interested in 
    
    G4String pfx_global = "/global/"; 

    fMessenger = new UserMessenger<RunParameters>(this); 

    //set which arm to use 
    fMessenger->AddCommand_string(
        pfx_global + "arm_mode",  //command name
        "arm_mode",               // parameter name
        &RunParameters::Set_ArmMode, // associated function in DetectorConstruction
        "both",                  // default value
        "both RHRS LHRS",             // possible valid inputs
        "Deterimes whether the LHRS or RHRS will be simulated. valid inputs are 'RHRS or LHRS'"
    ); 

    fMessenger->AddCommand_string(
        pfx_global + "outfile_path", 
        "outfile_path",
        &RunParameters::Set_OutfilePath,
        "test_output.root",
        "",
        "Determines path for output data file" 
    ); 

    fMessenger->AddCommand_string(
        pfx_global + "infile_path", 
        "infile_path",
        &RunParameters::Set_InfilePath,
        "none",
        "",
        "Determines path for input events file" 
    ); 

    //minimum momentum 
    fMessenger->AddCommand_doubleWithUnit(
        pfx_global + "momentum_low", 
        "momentum_low", 
        &RunParameters::Set_MomentumLow, 
        1104.*(1. - 0.06),
        "MeV", 
        "Minimum momentum of leptons to be saved to the output file"
    ); 

    //maximum momentum 
    fMessenger->AddCommand_doubleWithUnit(
        pfx_global + "momentum_high", 
        "momentum_high", 
        &RunParameters::Set_MomentumHigh, 
        1104.*(1. + 0.06),
        "MeV", 
        "Maximum momentum of leptons to be saved to the output file"
    ); 

    G4String pfx_generator = "/generator/";
    
    //minimum x-sieve to generate (mirrored between left / right arms)
    fMessenger->AddCommand_doubleWithUnit(
        pfx_generator + "sieve_x_low", 
        "sieve_x_low", 
        &RunParameters::Set_SieveXLow, 
        -85*mm,
        "mm",
        "Minimum x-val of target point on sieve-face"
    ); 

    //maximum x-sieve to generate (mirrored between left / right arms)
    fMessenger->AddCommand_doubleWithUnit(
        pfx_generator + "sieve_x_high", 
        "sieve_x_high", 
        &RunParameters::Set_SieveXHigh, 
        +85*mm,
        "mm",
        "Maximum x-val of target point on sieve-face"
    ); 
    
    //minimum y-sieve to generate (mirrored between left / right arms)
    fMessenger->AddCommand_doubleWithUnit(
        pfx_generator + "sieve_y_low", 
        "sieve_y_low", 
        &RunParameters::Set_SieveYLow, 
        -85*mm,
        "mm",
        "Minimum y-val of target point on sieve-face"
    ); 

    //maximum y-sieve to generate (mirrored between left / right arms)
    fMessenger->AddCommand_doubleWithUnit(
        pfx_generator + "sieve_y_high", 
        "sieve_y_high", 
        &RunParameters::Set_SieveYHigh, 
        +85*mm,
        "mm",
        "Maximum y-val of target point on sieve-face"
    ); 

    //vertical raster amplitude
    fMessenger->AddCommand_doubleWithUnit(
        pfx_generator + "vertical_raster_amplitude", 
        "vertical_raster_amplitude", 
        &RunParameters::Set_VerticalRasterAmplitude, 
        2*mm,
        "mm",
        "Vertical raster amplitude"
    ); 

    //horizontal raster amplitude
    fMessenger->AddCommand_doubleWithUnit(
        pfx_generator + "horizontal_raster_amplitude", 
        "horizontal_raster_amplitude", 
        &RunParameters::Set_HorizontalRasterAmplitude, 
        2*mm,
        "mm",
        "Horizontal raster amplitude"
    ); 

    //expected number of events to keep
    fMessenger->AddCommand_int(
        pfx_global + "expected_n_events_kept", 
        "expected_n_events_kept", 
        &RunParameters::Set_ExpectedNEventsKept, 
        0,
        "Horizontal raster amplitude"
    ); 

    //set which target mode to use  
    fMessenger->AddCommand_string(
        pfx_generator + "target_mode",  //command name
        "target_mode",               // parameter name
        &RunParameters::Set_TargetMode, // associated function in DetectorConstruction
        "production",                  // default value
        "production V1 V2 V3",             // possible valid inputs
        "Determines which target mode will be simulated (where to place react vertex)"
    ); 
}   
//______________________________________________________________________________
void RunParameters::Set_ArmMode(G4String arm_mode)
{
    const static std::map<G4String,ArmMode::EMode> val_map{
        {"both", ArmMode::kBoth},
        {"RHRS", ArmMode::kRHRS},
        {"LHRS", ArmMode::kLHRS}
    }; 

    const auto mode = val_map.find(arm_mode); 
    if (mode == val_map.end()) {

        G4Exception(
            Form("RunParameters::%s",__func__), 
            "Invalid arm mode passed", 
            G4ExceptionSeverity::RunMustBeAborted, 
            Form("Arm mode '%s' passed is invalid.",arm_mode.c_str())
        );
        return; 
    }
    f_arm_mode = mode->second; 
}
//______________________________________________________________________________
void RunParameters::Set_TargetMode(G4String arm_mode)
{
    const static std::map<G4String,TargetMode::Bit> val_map{
        {"production",  TargetMode::kProduction},
        {"V1",          TargetMode::kV1},
        {"V2",          TargetMode::kV2},
        {"V3",          TargetMode::kV3}
    }; 

    const auto mode = val_map.find(arm_mode); 
    if (mode == val_map.end()) {

        G4Exception(
            Form("RunParameters::%s",__func__), 
            "Invalid target mode passed", 
            G4ExceptionSeverity::RunMustBeAborted, 
            Form("Target mode '%s' passed is invalid.",arm_mode.c_str())
        );
        return; 
    }
    f_target_mode = mode->second; 
}
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________

