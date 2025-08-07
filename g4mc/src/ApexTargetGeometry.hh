#ifndef ApexTargetGeometry_HH
#define ApexTargetGeometry_HH

/////////////////////////////////////////////////////////////////////////////////////////
//
//  This is a quick and dirty way to make sure that any class who needs to know
//  specifics on the position / location of the sieve can have access to it. 
//
//  - Seth H. 8 Jul 25
//
/////////////////////////////////////////////////////////////////////////////////////////

#include "G4ThreeVector.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include <vector>

namespace ApexTargetGeometry {
  
  inline double Get_sieve_angle(bool _is_RHRS) { return ( _is_RHRS ? -5.369 : 5.369 ) * CLHEP::pi / 180.; }
  
  inline double Get_HRS_angle(bool _is_RHRS) { return ( _is_RHRS ? -12.50 : 12.50 ) * CLHEP::pi / 180.; }

  //units in mm. 
  inline G4ThreeVector Get_APEX_Target_center() { return G4ThreeVector( 0., 0., -1053.7952 ); }
    
  //units in mm. These are in Target coordinates (TCS), obtained by rotating hall
  // coordinates first by '-sieve_angle' about the y-axis, then by pi/2 about the
  // z-axis. 
  inline G4ThreeVector Get_sieve_pos(bool _is_RHRS) { return G4ThreeVector( 0., 0., 795.188 ); } 

  inline double Get_sieve_thickness() { return 12.7; }
  
  struct SieveHole {

    int row,col; 
    double x,y,radius_front,radius_back; 
    bool is_big; 
    
    //defining (overloading) this operator lets us use the std::find() function on a vector<SieveHole> 
    bool operator==(const SieveHole& rhs) { return ((row==rhs.row) && (col==rhs.col)); }
  }; 

  //constructs a container of SieveHole structs with accurate positions, row/column indices, and sizes. 
  // all units in mm for sieve hole positions and sizes. 
  std::vector<SieveHole> Construct_sieve_holes(bool _is_RHRS); 

}; 

#endif 
