#include "EventWriter.hh"
#include "G4Exception.hh"
#include "G4ExceptionSeverity.hh"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TParameter.h"

#include <stdexcept> 

//_____________________________________________________________________________
EventWriter::EventWriter() 
{
    /*noop*/ 
}
//_____________________________________________________________________________
void EventWriter::Initialize(ArmMode::EMode _mode, std::string path_outfile, long int max_events_estimate) 
{
    fArm_mode = _mode; 
    fData_singleArm.clear(); 
    fData_both.clear(); 
    if (fArm_mode == ArmMode::kBoth) { 
        fData_both.reserve(max_events_estimate); 
    } else {
        fData_singleArm.reserve(max_events_estimate); 
    }     

    fPath_outfile = path_outfile; 

    fIs_initialized = true; 
}; 
//_____________________________________________________________________________
void EventWriter::SaveEvent()
{
    switch (fArm_mode) {

        case ArmMode::kBoth : {
            fData_both.push_back({fTrackData_R, fTrackData_L}); 
            break; 
        } 
        case ArmMode::kRHRS : {
            fData_singleArm.push_back(fTrackData_R); 
            break; 
        } 
        case ArmMode::kLHRS : {
            fData_singleArm.push_back(fTrackData_L); 
            break; 
        } 
        default : {
            G4Exception(
                "EventWriter::SaveEvent",
                "EventWriter has invalid arm-mode",
                G4ExceptionSeverity::RunMustBeAborted,
                Form(
                    "'SaveEvent' was called, but EventWriter "
                    "static instance is not valid arm mode. EventWriter::GetArmMode() == %s", 
                ArmMode::Name(fArm_mode).c_str())
            ); 
            return; 
        }

    }
}
//_____________________________________________________________________________
void EventWriter::WriteOutputFile()
{
    auto file = new TFile(fPath_outfile.c_str(), "RECREATE");
    
    if (!file || !file->IsOpen()) {
        G4Exception(
            Form("EventWriter::%s",__func__), 
            "TFile unable to be opened",
            G4ExceptionSeverity::RunMustBeAborted, 
            Form("Unable create TFile under path '%s'", fPath_outfile.c_str()) 
        ); 
        return;
    }

    //
    auto Vec3_to_TVector3 = [](Vec3 v) { return TVector3(v.x, v.y, v.z); };
    //

    auto tree = new TTree("tracks_Q1", "Tracks propagated to Q1"); 

    if (fArm_mode == ArmMode::kBoth) {
        
        TrackData_t R_td, L_td; 

        TVector3 R_pos_vtx, R_mom_vtx, R_pos_Q1, R_mom_Q1, R_pos_sieve, R_mom_sieve; 
        TVector3 L_pos_vtx, L_mom_vtx, L_pos_Q1, L_mom_Q1, L_pos_sieve, L_mom_sieve; 

        //fill out the branches (generic, both-arm information)
        tree->Branch("event_id",         &R_td.event_id); 
        
        tree->Branch("target_code",      &R_td.target_code); 
        tree->Branch("target_num",       &R_td.target_num); 

        tree->Branch("invariant_mass",   &R_td.invariant_mass); 

        //arm-specific information
        //positron (RHRS)
        tree->Branch("R_track_id",       &R_td.track_id);
        
        tree->Branch("R_position_vtx",   &R_pos_vtx);
        tree->Branch("R_momentum_vtx",   &R_mom_vtx);
        
        tree->Branch("R_position_Q1",    &R_pos_Q1);
        tree->Branch("R_momentum_Q1",    &R_mom_Q1);
        
        tree->Branch("R_position_sieve", &R_pos_sieve);
        tree->Branch("R_momentum_sieve", &R_mom_sieve);

        //electron (LHRS)  
        tree->Branch("L_track_id",       &L_td.track_id);
        
        tree->Branch("L_position_vtx",   &L_pos_vtx);
        tree->Branch("L_momentum_vtx",   &L_mom_vtx);
        
        tree->Branch("L_position_Q1",    &L_pos_Q1);
        tree->Branch("L_momentum_Q1",    &L_mom_Q1);
        
        tree->Branch("L_position_sieve", &L_pos_sieve);
        tree->Branch("L_momentum_sieve", &L_mom_sieve);

        //loop over all data, write it to disk. 
        for (const auto& event : fData_both) {
            R_td = event.R; 
            L_td = event.L; 

            R_pos_Q1 = Vec3_to_TVector3(R_td.position_Q1);
            R_mom_Q1 = Vec3_to_TVector3(R_td.momentum_Q1);
            
            R_pos_vtx = Vec3_to_TVector3(R_td.position_vtx);
            R_mom_vtx = Vec3_to_TVector3(R_td.momentum_vtx);
            
            R_pos_sieve = Vec3_to_TVector3(R_td.position_sieve);
            R_mom_sieve = Vec3_to_TVector3(R_td.momentum_sieve);

            L_pos_Q1 = Vec3_to_TVector3(L_td.position_Q1);
            L_mom_Q1 = Vec3_to_TVector3(L_td.momentum_Q1);
            
            L_pos_vtx = Vec3_to_TVector3(L_td.position_vtx);
            L_mom_vtx = Vec3_to_TVector3(L_td.momentum_vtx);
            
            L_pos_sieve = Vec3_to_TVector3(L_td.position_sieve);
            L_mom_sieve = Vec3_to_TVector3(L_td.momentum_sieve);
            
            tree->Fill(); 
        }
        fData_both.clear(); 
    
    } else { 

        TrackData_t td; 
        TVector3 pos_vtx, mom_vtx, pos_Q1, mom_Q1, pos_sieve, mom_sieve; 

        //fill out the branches (generic, both-arm information)
        tree->Branch("event_id",         &td.event_id); 
        
        tree->Branch("target_code",      &td.target_code); 
        tree->Branch("target_num",       &td.target_num); 

        tree->Branch("invariant_mass",   &td.invariant_mass); 

        //arm-specific information
        //positron (RHRS)
        tree->Branch("track_id",       &td.track_id);
        
        tree->Branch("position_vtx",   &pos_vtx);
        tree->Branch("momentum_vtx",   &mom_vtx);
        
        tree->Branch("position_Q1",    &pos_Q1);
        tree->Branch("momentum_Q1",    &mom_Q1);
        
        tree->Branch("position_sieve", &pos_sieve);
        tree->Branch("momentum_sieve", &mom_sieve);
 
        //loop over all data, write it to disk. 
        for (const auto& event : fData_singleArm) {

            td = event; 

            pos_Q1 = Vec3_to_TVector3(td.position_Q1);
            mom_Q1 = Vec3_to_TVector3(td.momentum_Q1);
            
            pos_vtx = Vec3_to_TVector3(td.position_vtx);
            mom_vtx = Vec3_to_TVector3(td.momentum_vtx);
            
            pos_sieve = Vec3_to_TVector3(td.position_sieve);
            mom_sieve = Vec3_to_TVector3(td.momentum_sieve);
            
            tree->Fill(); 
        }
        fData_singleArm.clear(); 

        auto param_is_RHRS = new TParameter<bool>("is_RHRS", (fArm_mode==ArmMode::kRHRS)); 
        param_is_RHRS->Write(); 
    }

    tree->Write(); 
    file->Close(); 
}
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
