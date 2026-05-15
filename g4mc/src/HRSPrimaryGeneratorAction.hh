#ifndef HRSHRSPrimaryGeneratorAction_HH
#define HRSHRSPrimaryGeneratorAction_HH

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file HRSPrimaryGeneratorAction.hh
/// \brief Definition of the B1::HRSPrimaryGeneratorAction class

#include "G4VUserPrimaryGeneratorAction.hh"
#include "UserMessenger.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "ArmMode.hh"
#include "TargetMode.hh"

#include <random>

class G4ParticleGun;
class G4Event;
class G4Box;

/// The primary generator action class with particle gun.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued
/// in front of the phantom across 80% of the (X,Y) phantom size.

class HRSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  private:     
  bool f_is_initialized{false}; 
  public:
    HRSPrimaryGeneratorAction();
    ~HRSPrimaryGeneratorAction() override;

    // method from the base class
    void GeneratePrimaries(G4Event*) override;

  // method to access particle gun
  const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

  void Set_Verbose(G4int _x) { fVerbose=_x; }; 
  
  private:
    
    UserMessenger<HRSPrimaryGeneratorAction> *fMessenger; 

  G4int fVerbose{0}; 
  
    // center of the chosen target
    G4ThreeVector fTargetPosition; 

    G4double fBeamEnergy=2200.; 
    // vertical raster amplitude 
    G4double fRasterAmplitude_vertical; 

    G4RotationMatrix fSieveRotation; 

    // 'true' if RHRS, 'false' if LHRS
    G4bool f_is_RHRS; 

    ArmMode::EMode fArmMode; 

    inline bool RHRS_is_active() const { return fArmMode & ArmMode::kRHRS; }
    inline bool LHRS_is_active() const { return fArmMode & ArmMode::kLHRS; }
  
    TargetMode::Bit fTargetMode; 

    G4ThreeVector fGunPosition; 

    std::mt19937 fTwister; 
    std::uniform_int_distribution<int> fRandFoil{0,9};

    /// @return z-position of tungsten foil of target 
    inline double RandFoilZ(int& i_foil) {
      i_foil=fRandFoil(fTwister); 
      return -243.8 + ((double)i_foil)*55.0; 
    }

    G4double f_sieve_x_min, f_sieve_x_max, f_sieve_y_min, f_sieve_y_max; 
    G4double f_momentum_min, f_momentum_max; 
    G4double f_vertical_raster_amplitude, f_horizontal_raster_amplitude; 

    G4ParticleGun* fParticleGun = nullptr;  // pointer a to G4 gun class
    G4Box* fEnvelopeBox = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
