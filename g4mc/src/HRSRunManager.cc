#include "HRSRunManager.hh"
#include "G4String.hh"
#include "RunParameters.hh"

//initialize ptr to static singleton instance 
G4ThreadLocal HRSRunManager* HRSRunManager::fRunManager = nullptr; 

//_____________________________________________________________________________________________
HRSRunManager::HRSRunManager()
{
    
    //set the singleton instance to this object
    fRunManager = this; 

    fMessenger = new UserMessenger<HRSRunManager>(this); 

    G4String pfx_cmd = "/run/";

    //Set the input file to use (default is no input file) 
    fMessenger->AddCommand_string(
        pfx_cmd + "input_file",  //command name
        "input_file",            // parameter name
        &HRSRunManager::Set_InputFilePath, // associated function in DetectorConstruction
        "none",                  // default value
        "relative path to .root file of input events"
    ); 
}
//_____________________________________________________________________________________________
HRSRunManager::~HRSRunManager()
{
    if (fMessenger) delete fMessenger; 
}
//_____________________________________________________________________________________________
void HRSRunManager::Set_InputFilePath(G4String val)
{
    if (val=="none") {
        f_useInputFile=false;
        fInputFilePath="none"; 
    } else {
        f_useInputFile=true; 
        fInputFilePath=val; 
    }
}
//_____________________________________________________________________________________________
void HRSRunManager::BeamOn(G4int n_event, const char* macroFile, G4int n_select)
{
    //if we're using input events, then run **only** that many events. 
    if (UseInputFile()) {
        fEventReader = new EventReader(Get_InputFilePath()); 
        n_event = (G4int)fEventReader->GetNEvents(); //only read this many events from the event reader 
    }
    G4RunManager::BeamOn(n_event, macroFile, n_select); 
}
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________