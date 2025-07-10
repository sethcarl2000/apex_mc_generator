////////////////////////////////////////////////////////////////////////////////////////
//
//   This is a wrapper function which gives access to the relevant fortran subroutines.
//   Check simc_hrs.cc for detailed description of inputs/outputs.
//
//   The 'bool' reutrn value is 'true' if SIMC tells us the lepton has propagated thru
//   the HRS successfully, and 'false' if not. 
//
// !!! Definition of inputs to 'propagate_thru_hrs': !!!
//
//   x_tcs  -  Starting position of lepton in TARGET COORDINATE SYSTEM. Recall: TCS 
//   y_tcs     is the coordinate system in which x_tcs points down to the hall floor,
//   z_tcs     z points into the HRS opening, and y_hcs points to the left if you're
//             facing the HRS opening. (these arguments expect inputs in meters!)
//   dxdz_tcs
//   dydz_tcs  - Slopes in TCS. so: dxdz_tcs = (dx_tcs)/(dz_tcs), etc. (the Fortran
//               subroutine expects an input in slope-miliradians, so this script
//               will do the work of multiplying each slope by 1000. do not do this
//               beforehand!). 
//
//   dpp_tcs  - Defined as :  (p_track - p_HRS)/p_HRS, in which p_HRS is the central
//              momentum setting of the spectrometer (hard-coded in simc_hrs.hh), and
//              p_track is the momentum of the lepton before it enters the HRS. 
//     
//
//   Seth Hall 2 Jul 2025
//
////////////////////////////////////////////////////////////////////////////////////////

#include "simc_hrs.hh"
#include <cmath>
#include <memory>
#include <iostream>

using namespace std; 

extern "C" void mc_hrsr_(double *p_spec,
			 double *dpp,
			 double *x,
			 double *y,
			 double *z,
			 double *dxdz,
			 double *dydz,
			 double *x_fp,
			 double *dx_fp,
			 double *y_fp,
			 double *dy_fp,
			 double *m2,
			 bool *ms_flag,
			 bool *wcs_flag,
			 bool *decay_flag,
			 bool *skipto_Q1_flag,
			 double *fry,
			 bool *ok_spec,
			 double *pathlen); 

extern "C" void mc_hrsl_(double *p_spec,
			 double *dpp,
			 double *x,
			 double *y,
			 double *z,
			 double *dxdz,
			 double *dydz,
			 double *x_fp,
			 double *dx_fp,
			 double *y_fp,
			 double *dy_fp,
			 double *m2,
			 bool *ms_flag,
			 bool *wcs_flag,
			 bool *decay_flag,
			 bool *skipto_Q1_flag,
			 double *fry,
			 bool *ok_spec,
			 double *pathlen); 

bool SimcHandler::propagate_thru_hrs(SimcHandler::EArm arm,
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
				     double *path_length)
{

  double _p_spec(SimcHandler::hrs_momentum);
  double _dpp(dpp_tcs*100.);
  double _x(x_tcs*100.);
  double _y(y_tcs*100.);
  double _z(z_tcs*100.);
  double _dxdz(dxdz_tcs*1000.);
  double _dydz(dydz_tcs*1000.);
  double _m2(pow(0.510,2));
  bool   _ms_flag(false);
  bool   _wcs_flag(false);
  bool   _decay_flag(false);
  bool   _skipto_Q1_flag(true);
  double _fry(0.);
  bool   _ok_spec(false);
  
  if (arm==SimcHandler::kLeft) {
    mc_hrsl_(&_p_spec,
	     &_dpp,
	     &_x,
	     &_y,
	     &_z,
	     &_dxdz,
	     &_dydz,
	     x_fp,
	     dxdz_fp,
	     y_fp,
	     dydz_fp,
	     &_m2,
	     &_ms_flag,
	     &_wcs_flag,
	     &_decay_flag,
	     &_skipto_Q1_flag,
	     &_fry,
	     &_ok_spec,
	     path_length);
    cout << "done with s.r., result: " << (_ok_spec?"T":"F") << endl;
  } else                   {
    mc_hrsr_(&_p_spec,
	     &_dpp,
	     &_x,
	     &_y,
	     &_z,
	     &_dxdz,
	     &_dydz,
	     x_fp,
	     dxdz_fp,
	     y_fp,
	     dydz_fp,
	     &_m2,
	     &_ms_flag,
	     &_wcs_flag,
	     &_decay_flag,
	     &_skipto_Q1_flag,
	     &_fry,
	     &_ok_spec,
	     path_length);
  }
  
  return _ok_spec;
} 
