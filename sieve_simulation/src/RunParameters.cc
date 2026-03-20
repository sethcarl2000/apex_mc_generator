
#include "RunParameters.hh"
#include "G4String.hh"

namespace B1
{

//initialize ptr to the singleton instance 
RunParameters* RunParameters::fInstance = nullptr; 

//______________________________________________________________________________
RunParameters::RunParameters()
{
    //add the commands we're interested in 
    
    G4String cmd_prefix = "/global/"; 

    fMessenger = new UserMessenger<RunParameters>(this); 

    //set which arm to use 
    fMessenger->AddCommand_string(
        cmd_prefix + "active_arm",  //command name
        "active_arm",               // parameter name
        &RunParameters::SetArm, // associated function in DetectorConstruction
        "RHRS",                  // default value
        "RHRS LHRS",             // possible valid inputs
        "Deterimes whether the LHRS or RHRS will be simulated. valid inputs are 'RHRS or LHRS'"
    ); 

    fMessenger->AddCommand_string(
        cmd_prefix + "path_outfile", 
        "path_outfile",
        &RunParameters::SetPathOutfile,
        "output.root",
        "",
        "Determines path for output data file" 
    ); 

    //minimum momentum 
    fMessenger->AddCommand_doubleWithUnit(
        cmd_prefix + "min_momentum", 
        "min_momentum", 
        &RunParameters::SetMomentum_min, 
        1104.*(1. - 0.06),
        "MeV", 
        "Minimum momentum of leptons to be saved to the output file"
    ); 

    //maximum momentum 
    fMessenger->AddCommand_doubleWithUnit(
        cmd_prefix + "max_momentum", 
        "max_momentum", 
        &RunParameters::SetMomentum_max, 
        1104.*(1. + 0.06),
        "MeV", 
        "Maximum momentum of leptons to be saved to the output file"
    ); 
    

    G4String generator_prefix = "/generator/";

    //choose which target to use 
    fMessenger->AddCommand_string(
        generator_prefix + "target", 
        "target", 
        &RunParameters::SetTargetName, 
        "V2", 
        "V1 V2 V3", 
        "Deterimes which target to simulate"
    ); 

    //beam energy 
    fMessenger->AddCommand_doubleWithUnit(
        generator_prefix + "beam_energy", 
        "beam_energy", 
        &RunParameters::SetBeamEnergy, 
        2200.,
        "MeV", 
        "Beam energy"
    ); 
        
    fMessenger->AddCommand_string(
        generator_prefix + "mode", 
        "generator_mode",
        &RunParameters::SetGeneratorMode,
        "pair_production",
        "pair_production flat", 
        "Set the mode of e-/e+ generation for the particle gun"
    );

    fMessenger->AddCommand_doubleWithUnit(
        generator_prefix + "raster_amplitude",
        "raster_amplitude", 
        &RunParameters::SetVerticalRasterAmplitude,
        2.,
        "mm",
        "Set the (full) amplitude of the vertical raster"
    ); 
    
    fMessenger->AddCommand_string(
        cmd_prefix + "sieve_mode",
        "sieve_mode",
        &RunParameters::SetSieveMode,
        "all",
        "all small big wide_back",
        "Set which mode the sieve will be generated in"
    ); 

    fMessenger->AddCommand_doubleWithUnit(
        generator_prefix + "gunEnergy_min",
        "gunEnergy_min", 
        &RunParameters::SetGunEnergy_min,
        1000,
        "MeV",
        "Set the minimum particle-gun energy"
    ); 
    fMessenger->AddCommand_doubleWithUnit(
        generator_prefix + "gunEnergy_max",
        "gunEnergy_max", 
        &RunParameters::SetGunEnergy_max,
        2200,
        "MeV",
        "Set the maximum particle-gun energy"
    ); 

    fMessenger->AddCommand_doubleWithUnit(
        generator_prefix + "min_restMass",
        "min_restMass", 
        &RunParameters::SetMass_min,
        100.,
        "MeV",
        "Set the minimum of the decaying-particle rest-mass"
    ); 
    fMessenger->AddCommand_doubleWithUnit(
        generator_prefix + "max_restMass",
        "max_restMass", 
        &RunParameters::SetMass_max,
        400.,
        "MeV",
        "Set the maximum of the decaying-particle rest-mass"
    ); 
}   
//______________________________________________________________________________
void RunParameters::SetArm(G4String arm)
{
    if (arm == "RHRS") { f_is_RHRS = true; } else { f_is_RHRS = false; }
}
//______________________________________________________________________________
RunParameters* RunParameters::Instance()
{
    if (fInstance==nullptr) fInstance = new RunParameters(); 
    return fInstance; 
}
//______________________________________________________________________________
void RunParameters::SetGeneratorMode(G4String mode)
{ 
  if (mode == "pair_production")    { fGeneratorMode = kPairProduction; return; }
  if (mode == "flat")               { fGeneratorMode = kFlat; return; }
}
//______________________________________________________________________________
void RunParameters::SetSieveMode(G4String mode)
{
    if (mode == "all")          { fSieveMode = kAll; }
    if (mode == "small")        { fSieveMode = kSmall; }
    if (mode == "big")          { fSieveMode = kBig; }
    if (mode == "wide_back")    { fSieveMode = kWideBack; }
}
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________


}