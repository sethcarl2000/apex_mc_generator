#ifndef RunParameters_h 
#define RunParameters_h 1

#include "UserMessenger.hh"
#include "G4ParticleDefinition.hh"
#include "G4String.hh"

namespace B1 
{
    /// @brief enum to represent types of event generation
    enum EGeneratorMode { 
      kPairProduction=1, 
      kFlat
    }; 
    
    enum ESieveMode { kAll=0, kSmall, kBig, kWideBack }; 

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
    
    // minimum angle between the detected lepton and the beam for a track to be saved
    G4double fMinAngle;

    /// some parameters used by 'flat' generator ----------------------------------------------
    G4double fGunEnergy_max; 
    G4double fGunEnergy_min; 

    // the type of particle to generate (and throw at the target)
    G4ParticleDefinition* fGeneratedParticle; 

    // thickness of the target (unit - mm)
    G4double fTargetThickness; 


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
    void SetGeneratedParticle(G4String);
    void SetTargetThickness(G4double _x) { fTargetThickness=_x; }
    void SetGunEnergy_min(G4double _x) { fGunEnergy_min=_x; }
    void SetGunEnergy_max(G4double _x) { fGunEnergy_max=_x; }
    void SetMinAngle(G4double _x) { fMinAngle=_x; }
    //_________________________________________________________________________
    // 
    //   List of const methods to access read-only information 
    //

    /// @return 'true' if we're dealing with the RHRS, 'false' if LHRS 
    bool Is_RHRS() const { return f_is_RHRS; }
    
    /// @return path to output root file 
    G4String GetPathOutfile() const { return fPathOutfile; }

    /// @return minimum momentum of particles to save, MeV
    G4double GetMomentum_min() const { return fMomentum_min; }

    /// @return maximum momentum of particles to save, MeV
    G4double GetMomentum_max() const { return fMomentum_max; }
    
    /// @return minimum angle between outgoing particle and incident particle (to be saved)
    G4double GetMinAngle() const { return fMinAngle; }

    /// @return definition of particle to be thrown at target
    G4ParticleDefinition* GetParticleDefinition() const { return fGeneratedParticle; }
    
    /// @return thickness of target
    G4double GetTargetThickness() const { return fTargetThickness; }

    /// @return min particle gun energy
    G4double GetGunEnergy_min() const { return fGunEnergy_min; }

    /// @return max particle gun energy 
    G4double GetGunEnergy_max() const { return fGunEnergy_max; }

}; 

}

#endif 
