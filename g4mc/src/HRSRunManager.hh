#ifndef HRSStdHit_h_
#define HRSStdHit_h_

#include "G4RunManager.hh"
#include "UserMessenger.hh"
#include "G4String.hh"
#include "EventReader.hh"


// we need to overload this class so we can manage when the run generates its own events, and when it reads them from a file. 

class HRSRunManager : public G4RunManager {
private: 

    // Per-thread static instance of the run manager singleton.
    static G4ThreadLocal HRSRunManager* fRunManager;

    UserMessenger<HRSRunManager> *fMessenger; 

    bool f_useInputFile{false};
    G4String fInputFilePath{"none"}; 

    //event reader instance. run manager will own this. 
    EventReader* fEventReader{nullptr}; 

public: 
    
    //static access to singleton instance
    static HRSRunManager* GetRunManager() { return fRunManager; }

    //delete copy constructor & assignment operator
    HRSRunManager(const HRSRunManager&) = delete; 
    HRSRunManager& operator=(const HRSRunManager&) = delete; 

    //ordinary constructor & destructor 
    HRSRunManager(); 
    ~HRSRunManager();

    void Set_InputFilePath(G4String _x); 
    G4String Get_InputFilePath() const { return fInputFilePath; }

    /// @return 'true' if events are being read from an input file, 'false' otherwise 
    bool UseInputFile() const { return f_useInputFile; }

    //we're overriding this method to have it accept events from a ROOT file (if one is provided)
    void BeamOn(G4int n_event, const char* macroFile = nullptr, G4int n_select = -1) override; 

    /// @return ptr to event reader instance for this run. 
    EventReader* GetEventReader() { return fEventReader; } 
};


#endif