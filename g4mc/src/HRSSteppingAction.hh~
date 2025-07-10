// ********************************************************************
//
// $Id: HRSSteppingAction.hh,v2.0 2007/01/01 HRS Exp $
//
//..............................................................................

#ifndef HRSSteppingAction_h
#define HRSSteppingAction_h 1

#include "TFileHandler.hh"
#include "TrackData_t.hh"
#include "G4UserSteppingAction.hh"
#include "HRSSteppingActionMessenger.hh" 
#include "globals.hh"
#include <ostream>
#include <fstream>
#include <vector>
#include <memory>
#include "G4String.hh"
using namespace std;

//..............................................................................
class HRSSteppingAction : public G4UserSteppingAction
{
public:
  HRSSteppingAction();
  ~HRSSteppingAction();

  void UserSteppingAction(const G4Step*);
  
  inline void SetVerbose(G4int val) { verboseLevel = val; }
  inline G4int GetVerbose() const { return verboseLevel; }

  inline void EmptyPrintList() { vPrintList.clear(); }
  inline void Add2PrintList(G4String val) { vPrintList.push_back(val); }
  
  void SetOutFilePath(const std::string &_path);
  
  void SetOutTTreeName(const std::string &_name) {;} 

  void SetArm(HRSCoordinate::Arm _arm);
  HRSCoordinate::Arm GetArm() const;

  void Set_DoTXToutput(bool _val) { doTXToutput=_val; }
  bool Do_TXToutput() const { return doTXToutput; }

  void Set_killAtQ1(bool _val) { doKillAtQ1=_val; }
  bool Do_killAtQ1() const { return doKillAtQ1; }
  
private:
  void InitOutTxt();
  void DoPrint(const G4Step*);
  void PrintHead(const G4Step* theStep, ostream& pOut=G4cout);
  void PrintStep(const G4Step*, ostream& pOut=G4cout);
  void FillRootArray(const G4Step* theStep);
  
  bool write_trig;
  bool doTXToutput;
  bool doKillAtQ1;
  int evt_id, prevevt_id;
  std::ofstream output;
  std::ofstream output_vdc;
  std::ofstream output1;
  std::ofstream output_focpl;
  int i_st, iNoSecondary;

  //bool f_isRightArm; 
  
  //this is the 'TFileHandler' obj. which will let us
  //std::unique_ptr<TFileHandler> fOutFile; 
  //TFileHandler* fOutFile; 
  
  //std::unique_ptr<TrackData_t> fTrackData; 
  
  //std::unique_ptr<G4String> fOutFile_path;
  //G4String testStr;
  //std::unique_ptr<const char> fOutTTree_name; 
  //std::string fOutTTreeName; 

  HRSSteppingActionMessenger *messenger; 
  //std::unique_ptr<HRSSteppingActionMessenger> messenger;
  
  
  vector <G4String> vPrintList;
  bool CreateTxt;
  bool CreateRootNt;
  double mLHRSAngle;
  G4int mSeptumOn;
  double ang_fp;
  double P1x, P1y, P1z, E1v;
  double P2x, P2y, P2z, E2v;
  double px_tmp_mass, py_tmp_mass, pz_tmp_mass, e_tmp_mass;
  int el_n_mass, pos_n_mass;
  bool trig_gamma_conv;

  G4int verboseLevel;
  ofstream OutTxt;
  bool PrintHeadFlag;
  double z_vert_prev;

  int gam_evt_no;
  long int gam_evtNb;
  long int el_evtNb;
  long int pos_evtNb;
  double gam_z, gam_px, gam_py, gam_pz;
  double pos_z, pos_px, pos_py, pos_pz;
  double el_z, el_px, el_py, el_pz;
  bool  trig_pos_conv;
  bool  trig_el_conv;
  int  track_id_rad_tail_foil;
  double  px_rad_tail_foil, py_rad_tail_foil, pz_rad_tail_foil;
  int  track_id_sieve_back;
  double  px_sieve_back, py_sieve_back, pz_sieve_back;
  double  x_sieve_back, y_sieve_back, z_sieve_back;
};

//..............................................................................

#endif
