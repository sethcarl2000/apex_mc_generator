#include "EventReader.hh"

#include <G4Exception.hh>

// ROOT headers
#include <TString.h>
#include <ROOT/RDataFrame.hxx>

//stdlib headers
#include <stdexcept> 

//_____________________________________________________________________________________________
EventReader::EventReader(G4String path_input, G4String tree_name)
    : fPath_input{path_input}
{
    //try to open the input file 
    try {
    
        if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT();

        ROOT::RDataFrame df(tree_name.c_str(), path_input.c_str()); 
        
        fData = *df
            .Define("input_particle_info", [](
                double px, double py, double pz, 
                double rx, double ry, double rz, 
                char particle_type
            ) {
                return ParticleInput::Info_t{
                    .position_vtx = Vec3{ rx, ry, rz },
                    .momentum_vtx = Vec3{ px, py, pz },
                    .particle_type = particle_type
                }; 
            }, {
                "momentum_sieve_x", "momentum_sieve_y", "momentum_sieve_z", 
                "position_sieve_x", "position_sieve_y", "position_sieve_z"
                "particle_type"
            }).Take<ParticleInput::Info_t>("input_particle_info"); 
    
    } catch (const std::exception& e) {

        G4Exception(
            "EventReader::EventReader (constructor)", 
            "Exception caught trying to read input events",
            G4ExceptionSeverity::RunMustBeAborted, 
            Form("An exception was caught trying to read input events from file: '%s', what(): %s", e.what())
        );
        return; 
    }

    G4cout << "in <"<<__func__<<">: opened input file '"<< fPath_input <<"', with "<< fData.size() << " events read." << G4endl; 
}
//_____________________________________________________________________________________________
ParticleInput::Info_t EventReader::GetEvent(size_t i) const
{
    if (i >= fData.size() || i < 0) {
        G4Exception(
            "EventReader::GetEvent(int)", 
            "Invalid event number requested",
            G4ExceptionSeverity::EventMustBeAborted, 
            Form("Event number %i requested; valid range is [0,%i]", i, fData.size()-1)
        );
        return ParticleInput::Info_t{}; 
    }
    //if the event number requested was valid, return a copy of that particle 
    return fData[i];
}
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________

