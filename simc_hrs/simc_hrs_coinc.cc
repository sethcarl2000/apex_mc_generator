#include <stdio.h>
#include <iostream>
#include <stdio.h>

#include "simc_hrs.hh"
#include "fortran_functs/fortran_functs.hh"

using namespace std;

//extern "C" void my_subrt_(float *x,float *y,float *z, float *ans); 
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
			 double *fry,
			 double *pathlen,
			 bool *skipto_Q1_flag,
			 bool *ok_spec); 


//extern "C" void my_subrt_(float *x,float *y,float *z, float *ans); 
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
			 double *fry,
			 double *pathlen,
			 bool *skipto_Q1_flag,
			 bool *ok_spec); 

#include "HRSTrack_t.hh"

#include "ROOT/RDataFrame.hxx"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include "TStopwatch.h"
#include "TParameter.h"
//#include "RTypesCore.h"
#include <string>
#include <cstdio>
#include <memory>
#include <array>


//This is set in CMakeLists.txt now; to enable it, compile with -DENABLE_DEBUG=ON. 
#define DEBUG false

#define KEEP_FP_TRACKS_ONLY true

int main(int argc, char *argv[])
{
  //this is to disable implict MT... I don't trust Fortran to handle Implicit MT safely. 
  if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT(); 
  
  const double hrs_momentum = 1104.0; //in MeV/c
  
  const char* const here = "simc_hrs::main";
  
  //check if file was given. 
  if (argc<=1) {
    fprintf(stderr,
	    "<%s>: Error, no '.root' file was given as argument.\n\n"
	    "  usage: > simc_hrs 'path-to-my/file.root'\n",
	    here); 
    return 1; 
  }

  //use default outfile name, if none is given. 
  string path_outfile = "out_fp.root"; 
  string fp_tree_name = "tracks_fp"; 

  
  if (argc<=2) {
    printf("<%s>: Info: No output file name given. going with default: '%s'\n",
	   here, path_outfile.data());
  } else       {
    path_outfile = argv[2];
  }
    
  
  //now, attempt to open the input file
  string path_infile(argv[1]);
  
  auto infile = new TFile(path_infile.data(), "READ");
  
  //check if it was opened successfully
  if ( !infile ) {
    fprintf(stderr, "<%s>: Error: TFile '%s' could not be opened.\n",
	    here, path_infile.data()); 
    return 1;
  }
  if ( infile->IsZombie() ) {
    fprintf(stderr, "<%s>: Error: TFile '%s' is a zombie.\n",
	    here, path_infile.data()); 
    return 1;
  }
  
  
  //this is the name we expect for the input TTree
  const char* Q1_tree_name = "tracks_Q1";

  TTree *tree = (TTree*)infile->Get(Q1_tree_name); 

  //check if the tree could be found
  if (!tree) {
    fprintf(stderr, "<%s>: Error: Could not find TTree '%s' in input file '%s'.\n",
	    here, Q1_tree_name, path_infile.data());
    return 1; 
  }
  //get total number of events in the tree
  //long int n_events_read = tree->GetEntries(); 
  

  //now that we have checked that the TFile and desired TTree exist, we can open
  //the file and use RDataFrame.  
  
  infile->Close();
  delete infile;

  
  bool _skipto_Q1_flag(true); 
  
  //long int n_success(0); 
  
  bool _print_debug(false);
  /*#if DEBUG
    _print_debug=true;
  #endif*/ 
  
    
  //Now, open RDataFrame.
  ROOT::RDataFrame df(Q1_tree_name, path_infile.data());

  //________________________________________________________________________________________
  auto propagate_track_to_fp = [hrs_momentum,
				_print_debug,
				&_skipto_Q1_flag](TVector3 position_Q1,
						  TVector3 momentum_Q1, bool is_RHRS)
  {
    HRSTrack_t track;
      
    //track information at Q1 front window
    track.Q1.x = position_Q1.x();
    track.Q1.y = position_Q1.y();
    track.Q1.z = position_Q1.z();

    track.Q1.dxdz = momentum_Q1.x()/momentum_Q1.z();
    track.Q1.dydz = momentum_Q1.y()/momentum_Q1.z();
    
    track.Q1.dpp = (momentum_Q1.Mag() - hrs_momentum)/hrs_momentum;
      
    //this will be passed to the fortran subroutines.
#if DEBUG      
    printf("Q1-coords (% 4.1f,% 4.1f,% .5f,% .5f,% .5f)...",
	   track.Q1.x*100,
	   track.Q1.y*100,
	   track.Q1.dxdz,
	   track.Q1.dydz,
	   track.Q1.dpp);
#endif
      
    //now, we are ready to call 'mc_hrsl_/mc_hrsr_'. 
    double _p_spec(hrs_momentum);
    double _dpp(track.Q1.dpp * 100);
    double _x(track.Q1.x * 100.);
    double _y(track.Q1.y * 100.);
    double _z(172.050);
    double _dxdz(track.Q1.dxdz);
    double _dydz(track.Q1.dydz);
    double _x_tra(0.);
    double _dx_tra(0.);
    double _y_tra(0.);
    double _dy_tra(0.);
    double _m2(0.510*0.510);
    double _fry(0.);
    bool   _ok_spec(false);
    double _pathlen(0.);
      
    if (is_RHRS) { 
      mc_hrsr_(&_p_spec, &_dpp, &_x, &_y, &_z, &_dxdz, &_dydz,
	       &_x_tra, &_dx_tra, &_y_tra, &_dy_tra, &_m2, &_fry,
	       &_pathlen, &_skipto_Q1_flag, &_ok_spec);
    } else       {

      _y     *= -1.;
      _dydz  *= -1.; 

      //This should be the subroutine 'mc_hrsl_', but I found that there are some
      //acceptance issues with the simc_hrsl subroutine. It turns out that you get
      //good-looking results if you just feed the Left-arm data into 'mc_hrsr',
      //and use the fact that the spectrometers are mirrored (flip y & dy/dz). 
      mc_hrsr_(&_p_spec, &_dpp, &_x, &_y, &_z, &_dxdz, &_dydz,
	       &_x_tra, &_dx_tra, &_y_tra, &_dy_tra, &_m2, &_fry,
	       &_pathlen, &_skipto_Q1_flag, &_ok_spec);
	
      _y_tra  *= -1.;
      _dy_tra *= -1.; 
    }
      
#if DEBUG	
    cout << "print debug? " << (_print_debug?"True":"False") << endl; 
      
    cout << (_ok_spec?"pass":"fail");
    printf("pathlen:%-5.4f\n",_pathlen);
#endif

    //the fortran code outputs in cm, but we want the output in meters.
    _x_tra = _x_tra/100.;
    _y_tra = _y_tra/100.; 

      
    //Note: the fortran subroutine technically outputs coordinates in whats called
    // vdc 'transport' coordinates, hence the prefix 'tra'. to convert from
    // transport coordinates to focal-plane coordinates (fp), all you need to do is:
    //
    //    dxdz_fp = dxdz_tra - x_tra/6.
    // 
    track.fp.x    = _x_tra;
    track.fp.y    = _y_tra;
    track.fp.dxdz = _dx_tra - _x_tra/6.;
    track.fp.dydz = _dy_tra;

    track.flag_keep = _ok_spec; 

    track.path_length = _pathlen; 
      
    return track;
  }; 
  //________________________________________________________________________________________


  const auto n_events_total = (long int)*df.Count(); 
  
  
  auto df_output = df
    
    //.Range(0,1e7)
    
    .Define("R_track", [hrs_momentum, &propagate_track_to_fp]
	    (TVector3 position_vtx, TVector3 momentum_vtx,
	     TVector3 position_Q1,  TVector3 momentum_Q1)
    {
      //propagate track to fp
      auto track = propagate_track_to_fp(position_Q1, momentum_Q1, true); 

      //track information at vertex (target)
      track.tg.x = position_vtx.x();
      track.tg.y = position_vtx.y();
      track.tg.z = position_vtx.z();

      track.tg.dxdz = momentum_vtx.x()/momentum_vtx.z();
      track.tg.dydz = momentum_vtx.y()/momentum_vtx.z();
    
      track.tg.dpp = (momentum_vtx.Mag() - hrs_momentum)/hrs_momentum; 
      
      return track; 
      
    }, {"R_position_vtx", "R_momentum_vtx", "R_position_Q1", "R_momentum_Q1"})
  
    //Only keep an event if the track 'made it' to the focal plane.
    .Filter([](const HRSTrack_t &track){return track.flag_keep;}, {"R_track"})

    .Define("L_track", [hrs_momentum, &propagate_track_to_fp]
	    (TVector3 position_vtx, TVector3 momentum_vtx,
	     TVector3 position_Q1,  TVector3 momentum_Q1)
    {
      //propagate track to fp
      auto track = propagate_track_to_fp(position_Q1, momentum_Q1, false); 

      //track information at vertex (target)
      track.tg.x = position_vtx.x();
      track.tg.y = position_vtx.y();
      track.tg.z = position_vtx.z();

      track.tg.dxdz = momentum_vtx.x()/momentum_vtx.z();
      track.tg.dydz = momentum_vtx.y()/momentum_vtx.z();
    
      track.tg.dpp = (momentum_vtx.Mag() - hrs_momentum)/hrs_momentum; 
      
      return track; 
      
    }, {"L_position_vtx", "L_momentum_vtx", "L_position_Q1", "L_momentum_Q1"})
  
    //Only keep an event if the track 'made it' to the focal plane.
    .Filter([](const HRSTrack_t &track){return track.flag_keep;}, {"L_track"})
    
    //Deifine final output variables. This way, no one needs our class dictionary
    // to intepret the HRSTrack_t or HRSCoord_t structs
    .Define("R_x_sv",    [](const TVector3 &v){return v.x();},    {"R_position_sieve"})
    .Define("R_y_sv",    [](const TVector3 &v){return v.y();},    {"R_position_sieve"})
    .Define("R_dxdz_sv", [](const TVector3 &v){return v.x()/v.z();}, {"R_momentum_sieve"})
    .Define("R_dydz_sv", [](const TVector3 &v){return v.y()/v.z();}, {"R_momentum_sieve"})
    .Define("R_dpp_sv",  [hrs_momentum](const TVector3 &v) {return ( v.Mag() - hrs_momentum ) / hrs_momentum;}, {"R_momentum_sieve"})

    .Define("R_x_q1",    [](const TVector3 &v){return v.x();},    {"R_position_Q1"})
    .Define("R_y_q1",    [](const TVector3 &v){return v.y();},    {"R_position_Q1"})
    .Define("R_dxdz_q1", [](const TVector3 &v){return v.x()/v.z();}, {"R_momentum_Q1"})
    .Define("R_dydz_q1", [](const TVector3 &v){return v.y()/v.z();}, {"R_momentum_Q1"})
    .Define("R_dpp_q1",  [hrs_momentum](const TVector3 &v) {return ( v.Mag() - hrs_momentum ) / hrs_momentum;}, {"R_momentum_Q1"})
    
    .Define("R_x_fp",    [](const HRSTrack_t &track){return track.fp.x;}, {"R_track"})
    .Define("R_y_fp",    [](const HRSTrack_t &track){return track.fp.y;}, {"R_track"})
    .Define("R_dxdz_fp", [](const HRSTrack_t &track){return track.fp.dxdz;}, {"R_track"})
    .Define("R_dydz_fp", [](const HRSTrack_t &track){return track.fp.dydz;}, {"R_track"})
    .Define("R_is_at_fp",[](const HRSTrack_t &track){return track.flag_keep;},{"R_track"})
    .Define("R_path_length",[](const HRSTrack_t &track){return track.path_length;},{"R_track"})
    
    //Deifine final output variables. This way, no one needs our class dictionary
    // to intepret the HRSTrack_t or HRSCoord_t structs
    .Define("L_x_sv",    [](const TVector3 &v){return v.x();},    {"L_position_sieve"})
    .Define("L_y_sv",    [](const TVector3 &v){return v.y();},    {"L_position_sieve"})
    .Define("L_dxdz_sv", [](const TVector3 &v){return v.x()/v.z();}, {"L_momentum_sieve"})
    .Define("L_dydz_sv", [](const TVector3 &v){return v.y()/v.z();}, {"L_momentum_sieve"})
    .Define("L_dpp_sv",  [hrs_momentum](const TVector3 &v) {return ( v.Mag() - hrs_momentum ) / hrs_momentum;}, {"L_momentum_sieve"})

    .Define("L_x_q1",    [](const TVector3 &v){return v.x();},    {"L_position_Q1"})
    .Define("L_y_q1",    [](const TVector3 &v){return v.y();},    {"L_position_Q1"})
    .Define("L_dxdz_q1", [](const TVector3 &v){return v.x()/v.z();}, {"L_momentum_Q1"})
    .Define("L_dydz_q1", [](const TVector3 &v){return v.y()/v.z();}, {"L_momentum_Q1"})
    .Define("L_dpp_q1",  [hrs_momentum](const TVector3 &v) {return ( v.Mag() - hrs_momentum ) / hrs_momentum;}, {"L_momentum_Q1"})
    
    .Define("L_x_fp",    [](const HRSTrack_t &track){return track.fp.x;}, {"L_track"})
    .Define("L_y_fp",    [](const HRSTrack_t &track){return track.fp.y;}, {"L_track"})
    .Define("L_dxdz_fp", [](const HRSTrack_t &track){return track.fp.dxdz;}, {"L_track"})
    .Define("L_dydz_fp", [](const HRSTrack_t &track){return track.fp.dydz;}, {"L_track"})
    .Define("L_is_at_fp",[](const HRSTrack_t &track){return track.flag_keep;},{"L_track"})
    .Define("L_path_length",[](const HRSTrack_t &track){return track.path_length;},{"L_track"});
  
  printf("<%s>: Starting event loop...",here); cout << flush; 

  auto n_kept_rslt = df_output.Count(); 
  
  //measure the time it takes to go through all events
  TStopwatch timer; 

  //record all data to the final output tree
  df_output.Snapshot(fp_tree_name.data(), path_outfile.data(),
		     {"R_position_vtx",
		      "R_momentum_vtx",
		      "L_position_vtx",
		      "L_momentum_vtx",

		      "invariant_mass", 
			      
		      "R_x_sv",
		      "R_y_sv",
		      "R_dxdz_sv",
		      "R_dydz_sv",
		      "R_dpp_sv",

		      "R_x_q1",
		      "R_y_q1",
		      "R_dxdz_q1",
		      "R_dydz_q1",
		      "R_dpp_q1",

		      "R_x_fp",
		      "R_y_fp",
		      "R_dxdz_fp",
		      "R_dydz_fp",

		      "R_path_length",

		      "L_x_sv",
		      "L_y_sv",
		      "L_dxdz_sv",
		      "L_dydz_sv",
		      "L_dpp_sv",

		      "L_x_q1",
		      "L_y_q1",
		      "L_dxdz_q1",
		      "L_dydz_q1",
		      "L_dpp_q1",

		      "L_x_fp",
		      "L_y_fp",
		      "L_dxdz_fp",
		      "L_dydz_fp",

		      "L_path_length"});

  
  //done with event loop. 
  double time_elapsed = timer.RealTime(); 

  long int n_success     = (long int)*n_kept_rslt; 
  
  
  printf("\r<%s>: Starting event loop...Done.\n"
	 " --- Processed %li entries in %.3f seconds (%.3f ms/event)\n",
	 here, 
	 n_events_total,
	 time_elapsed,
	 1000.*time_elapsed/((double)n_events_total));
  
  printf(" --- Out of %li input events, %li were propagated (frac: %.4f)\n",
	 n_events_total,
	 n_success,
	 ((double)n_success)/((double)n_events_total));
  
  cout << flush; 
  
  return 0; 
}

  
