
#include "RunParameters.hh"
#include "G4ParticleTable.hh"
#include "G4String.hh"
#include "G4Exception.hh"
#include "G4ExceptionSeverity.hh"
#include "G4SystemOfUnits.hh"

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

    //path to output file 
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

    //minimum angle 
    fMessenger->AddCommand_doubleWithUnit(
        cmd_prefix + "min_angle", 
        "min_angle", 
        &RunParameters::SetMinAngle, 
        0.045,
        "rad", 
        "Minimum angle of scattered particles to save"
    ); 

    //thickness of the target 
    fMessenger->AddCommand_doubleWithUnit(
        "/detector/target_thickness", 
        "target_thickness", 
        &RunParameters::SetTargetThickness, 
        0.100,
        "mm", 
        "thickness of the target"
    ); 
        
    G4String generator_prefix = "/generator/";

    //name of particle to generate 
    fMessenger->AddCommand_string(
        generator_prefix + "particle_name", 
        "particle_name",
        &RunParameters::SetGeneratedParticle,
        "e-",
        "e- e+",
        "Determines which type of particle is thrown at the target" 
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
        generator_prefix + "gun_angle", 
        "gun_angle", 
        &RunParameters::SetGunAngle, 
        0*deg,
        "deg", 
        "Angle the particle gun makes with the normal"
    );

    //raster amplitude
    fMessenger->AddCommand_doubleWithUnit(
        generator_prefix + "raster_amplitude", 
        "raster_amplitude", 
        &RunParameters::SetRasterAmplitude, 
        20*mm,
        "mm", 
        "Amplitude of the 'Raster' (side-length of the square over which events will be generated)"
    ); 


}   
//______________________________________________________________________________
void RunParameters::SetGeneratedParticle(G4String particle_name)
{
    auto particle_table = G4ParticleTable::GetParticleTable(); 
    if (!particle_table) {
        G4Exception(__func__, "null ptr", G4ExceptionSeverity::RunMustBeAborted, "Ptr to Particle Table is null");
        fGeneratedParticle = nullptr; 
        return;  
    }
    fGeneratedParticle = particle_table->FindParticle(particle_name); 
    if (!fGeneratedParticle) {
        G4Exception(__func__, "particle not found", G4ExceptionSeverity::RunMustBeAborted, "Particle specified not found"); 
    }
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


}