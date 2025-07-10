#ifndef HRSTrack_t_HH
#define HRSTrack_t_HH

//#include "TObject.h"

////////////////////////////////////////////////////////////////////////////////////////////
//
// Very simple struct that we will stick in a vector to read/write to ROOT files
// Seth Hall 3 Jul 25
//
// Suggested use: all displacements in meters, and dpp defined as (p_trk - p_hrs)/p_hrs
//
////////////////////////////////////////////////////////////////////////////////////////////

struct HRSTrack_t {
public:

  struct HRSCoord_t {
    double x,y,z,dxdz,dydz,dpp;
  }; 
  
  //Coordinates at target, at Q1, and at fp. 
  HRSCoord_t tg,Q1,fp;

  bool flag_keep{false}; 


  double path_length; 
  //
    //double x_tg,y_tg,dxdz_tg,dydz_tg,dpp_tg;

  //double x_Q1,y_Q1,dxdz_Q1,dydz_Q1,dpp_Q1;

  //double x_fp,y_fp,dxdz_fp,dydz_fp;

}; 

#endif

