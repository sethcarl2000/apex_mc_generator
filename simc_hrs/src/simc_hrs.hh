#ifndef simc_hrs_HH
#define simc_hrs_HH

////////////////////////////////////////////////////////////////////////////////////////
//
//   This is a wrapper function which gives access to the relevant fortran subroutines.
//   Check simc_hrs.cc for detailed description of inputs/outputs
//
//   Seth Hall 2 Jul 2025
//
////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>

class SimcHandler {
public:

  SimcHandler() {}; 
  ~SimcHandler() {}; 
  
  enum EArm { kRight=0, kLeft }; 
  
  static constexpr double hrs_momentum = 1104.0;
  
  bool propagate_thru_hrs(SimcHandler::EArm arm,
			  double x_tcs,
			  double y_tcs,
			  double z_tcs,
			  double dxdz_tcs,
			  double dydz_tcs,
			  double dpp_tcs, 
			  double *x_fp,
			  double *y_fp,
			  double *dxdz_fp,
			  double *dydz_fp,
			  double *path_length);
}; 

#endif
