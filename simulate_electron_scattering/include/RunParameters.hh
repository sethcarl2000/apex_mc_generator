#ifndef RunParameters_h 
#define RunParameters_h 1

#include "UserMessenger.hh"
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
    G4String fTargetName; 
    G4double fBeamEnergy; 

    /// some parameters for the pair generator ------------------------------------------------
    // minimum rest mass for the decay particle
    G4double fMin_restMass;
    // maximum rest mass for the decay particle
    G4double fMax_restMass; 
    
    /// some parameters used by 'flat' generator ----------------------------------------------
    G4double fGunEnergy_max; 
    G4double fGunEnergy_min; 
    
    // (full) amplitude of the vertical raster 
    G4double fVerticalRasterAmplitude; 
    // center of target we're using
    G4ThreeVector fTargetPosition; 
    // the mode of paritlce generation 
    EGeneratorMode fGeneratorMode{kPairProduction}; 

    // the type of sieve generation (single-hole or all?)
    ESieveMode fSieveMode{kAll}; 

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
    void SetGeneratorMode(G4String mode); 
    void SetVerticalRasterAmplitude(G4double amplitude) { fVerticalRasterAmplitude=amplitude; }
    void SetMass_min(G4double _x) { fMin_restMass=_x; }
    void SetMass_max(G4double _x) { fMax_restMass=_x; }
    void SetSieveMode(G4String _v);

    void SetGunEnergy_min(G4double _x) { fGunEnergy_min=_x; }
    void SetGunEnergy_max(G4double _x) { fGunEnergy_max=_x; }
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

    /// @return the mode of particle generation
    EGeneratorMode GetGeneratorMode() const { return fGeneratorMode; }

    /// @return (full) amplitude of the vertical raster
    G4double GetRasterAmplitude_vertical() const { return fVerticalRasterAmplitude; }

    /// @return maximum mass of the pair generator
    G4double GetMass_max() const { return fMax_restMass; }

    /// @return minimum mass of the pair generator
    G4double GetMass_min() const { return fMin_restMass; }

    /// @return the sieve generator mode 
    ESieveMode GetSieveMode() const { return fSieveMode; }

    /// @return min particle gun energy
    G4double GetGunEnergy_min() const { return fGunEnergy_min; }

    /// @return max particle gun energy 
    G4double GetGunEnergy_max() const { return fGunEnergy_max; }
}; 

}

#endif 
