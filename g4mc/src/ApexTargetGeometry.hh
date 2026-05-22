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
  
  inline double Get_sieve_angle(bool _is_RHRS) {
    return ( _is_RHRS ? -5.372 : 5.366 ) * CLHEP::pi / 180.;
  }
  
  inline double Get_HRS_angle(bool _is_RHRS) {
    return ( _is_RHRS ? -12.50 : 12.50 ) * CLHEP::pi / 180.;
  }

  //units in mm. 
  inline G4ThreeVector Get_APEX_Target_center() {
    return G4ThreeVector( 0., 0., -1053.7952 );
  }
    
  //units in mm. These are in Target coordinates (TCS), obtained by rotating hall
  // coordinates first by '-sieve_angle' about the y-axis, then by pi/2 about the
  // z-axis. 
  inline G4ThreeVector Get_sieve_pos(bool _is_RHRS) {
    return G4ThreeVector(  _is_RHRS ?  -1.101 :  -1.301,
			   _is_RHRS ?  -3.885 :   6.672,
			   _is_RHRS ? 794.609 : 795.766 ); 
  } 

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

  /// @param v vector in HCS representing DISPLACEMENT from HCS origin 
  /// @param is_RHRS right(true) or left(false) arm
  /// @return vector in SCS representing displacement from SCS origin
  G4ThreeVector HCS_to_SCS(const G4ThreeVector& v, bool is_RHRS); 

  /// @param v vector in HCS representing DISPLACEMENT from SCS origin 
  /// @param is_RHRS right(true) or left(false) arm
  /// @return vector in SCS representing displacement from HCS origin
  G4ThreeVector SCS_to_HCS(const G4ThreeVector& v, bool is_RHRS); 

  /// @brief Project a position & displacement (given in HCS) onto the face of the sieve
  /// @param is_RHRS right(true) or left(false) HRS
  /// @param R displacement from center of HCS (center of APEX scattering chamber)
  /// @param S vector representing direction in HCS
  /// @param x x_sv
  /// @param y y_sv
  /// @param dxdz dx/dz_sv
  /// @param dydz dy/dz_sv
  void Project_HCS_onto_sieve(bool is_RHRS, const G4ThreeVector& R, const G4ThreeVector& S, double &x, double &y, double& dxdz, double& dydz); 

  /// @brief Given a position on the face of the sieve, return that point's position in HCS 
  /// @param is_RHRS right / left arm 
  /// @param x_sv x-pos on sieve face 
  /// @param y_sv y-pos on sieve face 
  /// @return Displacement from this point on the face of the sieve relative to the APEX scattering chamber center. 
  inline G4ThreeVector Get_sieve_intercept_HCS(bool is_RHRS, double x_sv, double y_sv) {
    return SCS_to_HCS(G4ThreeVector(x_sv, y_sv, 0.), is_RHRS); 
  };

}; 

#endif 
