#ifndef TFILEHANDLER_HH
#define TFILEHANDLER_HH

//////////////////////////////////////////////////////////////////////////////////////
//  
//  This is a placeholder class meant to be able to store the tracks 
//  in a ROOT file, rather than ASCII output.
//
//  Most data is collected in HRSSteppingAction.cc.
//  
//////////////////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TObject.h" 
#include "TFile.h"
#include "TTree.h"
#include <memory>
#include "TParameter.h"

#include "TrackData_t.hh"

class TFileHandler : public TObject {
public: 
  
  TFileHandler(const char *const _path,
	       bool is_RHRS, 
	       const char *const _treeName="tracks_Q1");

  ~TFileHandler();  
  
  void CloseFile();

  void ClearTrack(bool write_track=false); 
  
  TrackData_t* Get_TData_p() { return &fTrackData_p; }  
  TrackData_t* Get_TData_e() { return &fTrackData_e; }
  
  void WriteTrack(); 
  
  enum ETrackStatus { kNone=0, kDead, kAlive, kQ1 };

  ETrackStatus GetStatus(bool is_RHRS) const { return is_RHRS ? fStatus_p : fStatus_e; }

  void SetStatus(bool is_RHRS, ETrackStatus stat) {
    if (is_RHRS) { fStatus_p=stat; } else { fStatus_e=stat; }
  }
  
private:
  
  TrackData_t fTrackData_p; 
  TrackData_t fTrackData_e; 
  
  std::unique_ptr<TFile> fFile; 
  
  TTree *fTree; 

  TParameter<bool> *fParam_isRHRS; 

  ETrackStatus fStatus_p, fStatus_e; 
  
  ClassDef(TFileHandler,1); 
};
#endif
