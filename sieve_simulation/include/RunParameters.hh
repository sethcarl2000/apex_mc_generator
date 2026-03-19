#ifndef RunParameters_h 
#define RunParameters_h 1

#include "UserMessenger.hh"
#include "G4String.hh"

namespace B1 
{
    

/// @brief Singleton class which stores global run parameters, accesible by any class
class RunParameters {
private: 
    
    //the singleton instance 
    static RunParameters *fInstance; 

    /// The associated messenger
    UserMessenger<RunParameters> *fMessenger{nullptr}; 

    ///now, list all parameters of interest 
    G4bool f_is_RHRS; 
    G4String fPathOutfile;
    G4double fMomentum_min, fMomentum_max; 
    G4String fTargetName; 
    G4double fBeamEnergy; 

public: 

    RunParameters(); 
    ~RunParameters() {}; 

    /// @return ptr to singleton instance of global RunParameters object 
    static RunParameters* Instance(); 

    
    //_________________________________________________________________________
    // 
    //   List setters to set run information 
    //
    void SetArm(G4String);
    void SetPathOutfile(G4String _val)  { fPathOutfile=_val; } 
    void SetMomentum_min(G4double _val) { fMomentum_min=_val; }
    void SetMomentum_max(G4double _val) { fMomentum_max=_val; }
    void SetTargetName(G4String _val)   { fTargetName=_val; }
    void SetBeamEnergy(G4double _val)   { fBeamEnergy=_val; }

    //_________________________________________________________________________
    // 
    //   List of const methods to access read-only information 
    //

    /// @return 'true' if we're dealing with the RHRS, 'false' if LHRS 
    bool Is_RHRS() const { return f_is_RHRS; }
    
    /// @return path to output root file 
    G4String GetPathOutfile() const { return fPathOutfile; }

    /// @return minimum momentum, MeV
    G4double GetMomentum_min() const { return fMomentum_min; }

    /// @return maximum momentum, MeV
    G4double GetMomentum_max() const { return fMomentum_max; }

    /// @return name of target
    G4String GetTargetName() const { return fTargetName; }

    /// @return beam energy (MeV)
    G4double GetBeamEnergy() const { return fBeamEnergy; }
}; 

}

#endif 
