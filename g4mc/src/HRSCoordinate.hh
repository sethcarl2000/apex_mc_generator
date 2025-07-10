#ifndef HRSCoordinate_HH
#define HRSCoordinate_HH

////////////////////////////////////////////////////////////////////////////////////////
//
//  Very simple struct meant to keep a pair of G4ThreeVector's, one for a position,
//  and one for a momentum
//
//
//
////////////////////////////////////////////////////////////////////////////////////////

#include "G4ThreeVector.hh"

struct HRSCoordinate {
  
public:
  
  G4ThreeVector position { G4ThreeVector(0,0,0) };
  G4ThreeVector momentum { G4ThreeVector(0,0,0) };
  
  enum Coord { kHCS,   kTCS }; 
  enum Arm   { kNone, kRight, kLeft }; 
  
  Coord coord { kHCS };
  Arm arm { kNone }; 
}; 

#endif 
