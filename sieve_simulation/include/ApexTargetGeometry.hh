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
#include "G4SystemOfUnits.hh"

namespace ApexTargetGeometry {
  
  inline double Get_sieve_angle(bool _is_RHRS) {
    return ( _is_RHRS ? -5.372*deg : 5.366*deg );
  }
  
  inline double Get_HRS_angle(bool _is_RHRS) {
    return ( _is_RHRS ? -12.50*deg : 12.50*deg );
  }

  //units in mm. 
  inline G4ThreeVector Get_APEX_Target_center() {
    return G4ThreeVector( 0.*mm, 0.*mm, -1053.7952*mm );
  }
    
  //units in mm. These are in Target coordinates (TCS), obtained by rotating hall
  // coordinates first by '-sieve_angle' about the y-axis, then by pi/2 about the
  // z-axis. 
  inline G4ThreeVector Get_sieve_pos(bool _is_RHRS) {
    return G4ThreeVector(  
        _is_RHRS ?  -1.101*mm :  -1.301*mm,
		_is_RHRS ?  -3.885*mm :  +6.672*mm,
		_is_RHRS ? 794.609*mm : 795.766*mm 
    ); 
  } 
  
}; 

#endif 