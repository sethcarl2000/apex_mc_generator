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

struct TrackData_t {

public:

  enum EParticleType { kNone=0, kElectron, kPositron };
  
  EParticleType particle_type { EParticleType::kNone }; 
  
  HRSCoordinate::Arm arm { HRSCoordinate::Arm::kNone };
  
  int event_id {-1};
  int track_id {-1};
  
  //position & momentum at sieve plane
  //HRSCoordinate coord_sieve_HCS {}; 
  
  //position & momentum at Q1_front (always HCS)
  //HRSCoordinate coord_Q1 {};  

  //NOTE: these are both in 'rotated' hall coorindates;
  //for each arm, they are rotated as if the spectrometer they're entering is at 0-deg.
  TVector3 position_Q1;
  TVector3 momentum_Q1;

  //These are just in standard hall coordinates. 
  TVector3 position_vtx;
  TVector3 momentum_vtx;

  //These are in TARGET COORDINATES, projected to the sieve plane. 
  TVector3 position_sieve;
  TVector3 momentum_sieve; 
  
  //position & momentum at vertex
  //HRSCoordinate coord_vertex {}; 
};

#endif
