#ifndef EventReader_HH
#define EventReader_HH

#include "TrackData_t.hh"
#include "ArmMode.hh"
#include "UserMessenger.hh"
#include <string> 
#include <vector>

namespace ParticleInput {

    constexpr char kElectron = 'e'; 
    constexpr char kPositron = 'p'; 

    struct Info_t {
        Vec3 position_vtx, momentum_vtx; 
        char particle_type=kElectron; 
    }; 
};

//class to read input particles from an input file
class EventReader {
private: 
    
    //list of input events
    std::vector<ParticleInput::Info_t> fData; 
    
    std::string fPath_input; 

    ArmMode::EMode fArm_mode{ArmMode::kNone}; 

public:
    
    EventReader(G4String path_input, G4String tree_name="tracks_sieve"); 
    ~EventReader() {};

    //attempt to get access to event 'i' 
    ParticleInput::Info_t GetEvent(size_t i) const; 

    size_t GetNEvents() const { return fData.size(); }; 
};

static_assert(std::is_trivially_copy_constructible_v<ParticleInput::Info_t>, "ParticleInput::Info_t is not trivially copy constructable"); 
static_assert(std::is_trivially_move_constructible_v<ParticleInput::Info_t>, "ParticleInput::Info_t is not trivially move constructable"); 


#endif