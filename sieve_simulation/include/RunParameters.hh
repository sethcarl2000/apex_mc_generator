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


    //_________________________________________________________________________
    // 
    //   List of const methods to access read-only information 
    //

    /// @return 'true' if we're dealing with the RHRS, 'false' if LHRS 
    bool Is_RHRS() const { return f_is_RHRS; }
    
}; 

}

#endif 
