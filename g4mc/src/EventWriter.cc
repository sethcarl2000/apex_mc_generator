#include "EventWriter.hh"
#include "G4Exception.hh"
#include "G4ExceptionSeverity.hh"
#include "TString.h"

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
    /*noop*/
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
