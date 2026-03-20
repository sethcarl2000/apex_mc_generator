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

  G4String target_name = run_params->GetTargetName(); 
  
  //set the electron beam generation spot to be a little upstream
  fTargetPosition  
    = ApexTargetGeometry::GetTargetPosition(target_name); 

  fRasterAmplitude_vertical = run_params->GetRasterAmplitude_vertical(); 

  fBeamEnergy = run_params->GetBeamEnergy(); 

  fGeneratorMode = run_params->GetGeneratorMode(); 

  fMin_restMass = run_params->GetMass_min();
  fMax_restMass = run_params->GetMass_max();

  f_is_RHRS = run_params->Is_RHRS(); 

  if (fGeneratorMode == kPairProduction) {
    if (run_params->Is_RHRS()) {
      fParticleGun -> SetParticleDefinition(G4Positron::Positron()); 
    } else {
      fParticleGun -> SetParticleDefinition(G4Electron::Electron()); 
    }
  }    
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
  using namespace std; 
  
  // beam energy 
  G4double beam_E = fBeamEnergy; 

  // pick a random decay invariant mass
  G4double mass = fMin_restMass  +  (fMax_restMass - fMin_restMass)*G4UniformRand(); 

  // cosine between decay particle direction and beam (in the restframe of the particle)
  G4double cosine_restframe = (1 - 2*G4UniformRand());  

  // azimuthal angle of lepton around beam (defined so that phi=0 means it's in the horizontal plane)
  G4double phi = kPhi_max * (1 - 2*G4UniformRand()); 

  G4double sine_restframe = sqrt( 1. - cosine_restframe*cosine_restframe ); 

  // we're pretending the lepton is massless. 
  auto p_lepton = G4ThreeVector( 
    cos(phi) * sine_restframe * (f_is_RHRS ? -1 : +1), 
    sin(phi) * sine_restframe, 
    cosine_restframe
  );
  p_lepton *= (mass/2);

  // gamma factor of decaying particle in lab frame
  G4double gamma = beam_E/mass; 

  // boost the lepton to the lab frame
  p_lepton[2] = gamma*( p_lepton[2] + (mass/2)*sqrt(1. - 1./(gamma*gamma)) );  

  fParticleGun->SetParticleEnergy(p_lepton.mag()); 
  fParticleGun->SetParticleMomentumDirection(p_lepton); 

  // now, pick a random reaction vertex at the target
  fParticleGun->SetParticlePosition(
    G4ThreeVector(
      fTargetPosition.x(), 
      fTargetPosition.y() + (1 - 2*G4UniformRand())*fRasterAmplitude_vertical,
      fTargetPosition.z() 
    )
  );
      
  fParticleGun->GeneratePrimaryVertex(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
