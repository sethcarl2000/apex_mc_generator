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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"
#include "RunParameters.hh"
#include "ApexTargetGeometry.hh"
#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "Randomize.hh"

#include <cmath> 

namespace
{ 
  // maximum value of phi (azimuth around beam) from the horizontal of a generated lepton 
  constexpr G4double kPhi_max = 90.*deg; 

  // minimum angle between the beam and a generated lepton 
  constexpr G4double kTheta_min = 0.045*rad; 
} // namespace


namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  fMessenger = new UserMessenger<PrimaryGeneratorAction>(this); 

  G4String command_prefix = "/generator/"; 


  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  auto run_params = RunParameters::Instance(); 

  // default particle kinematic
  //fParticleGun->SetParticleDefinition(G4Electron::Electron());
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  //fParticleGun->SetParticleEnergy(run_params->GetBeamEnergy());
  
  f_is_RHRS = run_params->Is_RHRS(); 

  fParticleGun->SetParticleDefinition(run_params->GetParticleDefinition()); 
  
  fBeamEnergy_min = run_params->GetGunEnergy_min(); 
  fBeamEnergy_max = run_params->GetGunEnergy_max(); 

  fGunAngle = run_params->GetGunAngle(); 
  fRasterAmplitude = run_params->GetRasterAmplitude(); 

  //our 'pivot point' we rotate the gun around is the back of the target
  G4double dist_to_pivot = 1*cm + 3*mm + run_params->GetTargetThickness(); 
  
  //set the angle of the particle gun
  fParticleGun->SetParticleMomentumDirection(
    G4ThreeVector(std::sin(fGunAngle), 0, std::cos(fGunAngle))
  ); 
  
  fGunStartingPos = G4ThreeVector(
    -dist_to_pivot*std::tan(fGunAngle),
    0, 
    -run_params->GetTargetThickness()/2 -1*cm 
  ); 

  fParticleGun->SetParticlePosition(fGunStartingPos);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{ 
  //choose a random incident particle energy, in a random user-defined range
  fParticleGun->SetParticleEnergy(
    (fBeamEnergy_max - fBeamEnergy_min)*G4UniformRand() + fBeamEnergy_min
  );

  //choose a random starting point to launch the electron from
  fParticleGun->SetParticlePosition(G4ThreeVector(
    fGunStartingPos.x() + (2. - 1.)*G4UniformRand()*fRasterAmplitude, 
    fGunStartingPos.y() + (2. - 1.)*G4UniformRand()*fRasterAmplitude, 
    fGunStartingPos.z()
  )); 

  fParticleGun->GeneratePrimaryVertex(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
