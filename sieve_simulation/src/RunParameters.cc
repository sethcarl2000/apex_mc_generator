
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

    fMessenger->AddCommand_doubleWithUnit(
        cmd_prefix + "min_momentum", 
        "min_momentum", 
        &RunParameters::SetMomentum_min, 
        1104.*(1. - 0.06),
        "MeV", 
        "Minimum momentum of leptons to be saved to the output file"
    ); 

    fMessenger->AddCommand_doubleWithUnit(
        cmd_prefix + "max_momentum", 
        "max_momentum", 
        &RunParameters::SetMomentum_max, 
        1104.*(1. + 0.06),
        "MeV", 
        "Maximum momentum of leptons to be saved to the output file"
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
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________


}