#ifndef TrackData_t_HH
#define TrackData_t_HH

//////////////////////////////////////////////////////////////////////////////////////
//  
//   This is a pretty straightforward struct which will be created at the start
//   of each event within HRSSteppingAction::UserSteppingAction. It will record
//   all data needed to create the output root file. 
//  
//////////////////////////////////////////////////////////////////////////////////////

#include "G4ThreeVector.hh"
#include "HRSCoordinate.hh"
#include "TROOT.h"
#include "TVector3.h"
#include "Vec3.hh"
#include "TargetMode.hh"
#include <array> 

struct TrackData_t {

  enum EParticleType { kNone=0, kElectron, kPositron };
  
  EParticleType particle_type { EParticleType::kNone }; 

  enum EStatus { kDead=0, kAlive, kQ1 };
  EStatus status{EStatus::kDead};
  
  HRSCoordinate::Arm arm { HRSCoordinate::Arm::kNone };
  
  int event_id {-1};
  int track_id {-1};

  //target code. 
  // 0 - no target specified. 
  // W - production foils 
  // V - vertical wires
  // H - horizontal wires 
  // O - optics foils
  char target_code='0';
  
  // index of target. ex;  'V1' would have code 'V' and num '1'. 
  int target_num;



  //position & momentum at sieve plane
  //HRSCoordinate coord_sieve_HCS {}; 
  
  //position & momentum at Q1_front (always HCS)
  //HRSCoordinate coord_Q1 {};  

  //NOTE: these are both in 'rotated' hall coorindates;
  //for each arm, they are rotated as if the spectrometer they're entering is at 0-deg.
  //TVector3 position_Q1;
  //TVector3 momentum_Q1;
  Vec3 position_Q1, momentum_Q1; 


  //These are just in standard hall coordinates. 
  //TVector3 position_vtx;
  //TVector3 momentum_vtx;
  Vec3 position_vtx, momentum_vtx; 

  //invariant mass of event that generated this lepton
  double invariant_mass;
  
  //These are in TARGET COORDINATES, projected to the sieve plane. 
  Vec3 position_sieve, momentum_sieve; 
  //TVector3 position_sieve;
  //TVector3 momentum_sieve; 
  
  //position & momentum at vertex
  //HRSCoordinate coord_vertex {}; 
};

static_assert(std::is_trivially_copy_constructible_v<TrackData_t>, "TrackData_t is not trivially copy constructable"); 
static_assert(std::is_trivially_move_constructible_v<TrackData_t>, "TrackData_t is not trivially move constructable"); 

#endif
