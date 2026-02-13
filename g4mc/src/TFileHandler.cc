//////////////////////////////////////////////////////////////////////////////////////
//  
//  This is a placeholder class meant to be able to store the tracks 
//  in a ROOT file, rather than ASCII output.
//
//  Most data is collected in HRSSteppingAction.cc.
//  
//////////////////////////////////////////////////////////////////////////////////////

#include "TFileHandler.hh"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "G4ThreeVector.hh"
//#inlcude "TParameter.h"
#include <memory>

using namespace std; 

//____________________________________________________________________________________
TFileHandler::TFileHandler(const char *const _path, bool is_RHRS, const char *const _treeName)
  : fTrackData_p(TrackData_t{}), fTrackData_e(TrackData_t{})
{
  fFile = unique_ptr<TFile>(new TFile(_path,"RECREATE"));
  //fTree = unique_ptr<TTree>(new TTree(_treeName,"Output tree - Q1/vtx data"));
  fTree = new TTree(_treeName,"Output tree - Q1/vtx data");
  fTree->SetMakeClass(1);
  
  //fill out the branches:
  fTree->Branch("event_id",         &fTrackData_p.event_id); 

  
  //positron (RHRS)
  fTree->Branch("R_track_id",       &fTrackData_p.track_id);
  
  fTree->Branch("R_position_vtx",   &fTrackData_p.position_vtx);
  fTree->Branch("R_momentum_vtx",   &fTrackData_p.momentum_vtx);
  
  fTree->Branch("R_position_Q1",    &fTrackData_p.position_Q1);
  fTree->Branch("R_momentum_Q1",    &fTrackData_p.momentum_Q1);
  
  fTree->Branch("R_position_sieve", &fTrackData_p.position_sieve);
  fTree->Branch("R_momentum_sieve", &fTrackData_p.momentum_sieve);

  
  //electron (LHRS)  
  fTree->Branch("L_track_id",       &fTrackData_e.track_id);
  
  fTree->Branch("L_position_vtx",   &fTrackData_e.position_vtx);
  fTree->Branch("L_momentum_vtx",   &fTrackData_e.momentum_vtx);
  
  fTree->Branch("L_position_Q1",    &fTrackData_e.position_Q1);
  fTree->Branch("L_momentum_Q1",    &fTrackData_e.momentum_Q1);
  
  fTree->Branch("L_position_sieve", &fTrackData_e.position_sieve);
  fTree->Branch("L_momentum_sieve", &fTrackData_e.momentum_sieve);

  
  fTree->Branch("invariant_mass", &fTrackData_p.invariant_mass); 
  
  fParam_isRHRS = new TParameter<bool>("is_RHRS", is_RHRS);
  fParam_isRHRS->Write(); 
}
//____________________________________________________________________________________
TFileHandler::~TFileHandler()
{
  if (fFile) { if (fFile && fFile->IsOpen()) CloseFile(); }
  //delete fTree;
  //delete fFile; 
}
//____________________________________________________________________________________
void TFileHandler::WriteTrack()
{
  //this class will handle filling of the TTree
  fTree->Fill(); 
  return;
}
//____________________________________________________________________________________
void TFileHandler::CloseFile()
{
  if (fFile && fFile->IsOpen()) {
    if (!fFile->IsZombie()) fTree->Write(); 
    delete fTree; 
    delete fParam_isRHRS; 
    fFile->Close();
  }
}
//____________________________________________________________________________________
void TFileHandler::ClearTrack(bool write_track)
{
  //fTrackData = TrackData_t{}; 
  

}
//____________________________________________________________________________________
//____________________________________________________________________________________
//____________________________________________________________________________________
//____________________________________________________________________________________
//____________________________________________________________________________________
ClassImp(TFileHandler);
