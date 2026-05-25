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
#include "EnergyAngleGeneratorFactory.hh"

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
#include "G4Exception.hh"

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
  
  switch (run_params->GetGeneratorMode()) {

    case kElastic           : fEnergyAngleGenerator = EnergyAngleGeneratorFactory::Elastic(); break; 
    case kTrident_electron  : fEnergyAngleGenerator = EnergyAngleGeneratorFactory::Trident_Electron(); break; 
    case kTrident_positron  : fEnergyAngleGenerator = EnergyAngleGeneratorFactory::Trident_Positron(); break; 
    case kBHPhotoproduction : fEnergyAngleGenerator = EnergyAngleGeneratorFactory::BetheHeitlerPhotoproduction(); break; 
    default : {
      G4Exception(
        "PrimaryGeneratorAction::PrimaryGeneratorAction", 
        "Unsupported particle generator type",
        G4ExceptionSeverity::RunMustBeAborted, 
        "Particle type gotten from RunParameters is unsupported."
      );
    }
  }

  fParticleGun->SetParticleDefinition(fEnergyAngleGenerator->GetDefinition()); 

  //set the electron beam generation spot to be a little upstream
  fTargetPosition  
    = ApexTargetGeometry::GetTargetPosition(target_name); 

  fRasterAmplitude_vertical = run_params->GetRasterAmplitude_vertical(); 

  fBeamEnergy = run_params->GetBeamEnergy(); 

  fGeneratorMode = run_params->GetGeneratorMode(); 

  fMin_restMass = run_params->GetMass_min();
  fMax_restMass = run_params->GetMass_max();

  f_is_RHRS = run_params->Is_RHRS(); 

  fSieveRotation = 
    CLHEP::HepRotationY(ApexTargetGeometry::Get_sieve_angle(f_is_RHRS)) * CLHEP::HepRotationZ(-CLHEP::pi/2); 
  

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
  
  auto run_params = RunParameters::Instance(); 

  //pick a reaction vertex
  G4ThreeVector vertex(
    fTargetPosition.x(), 
    fTargetPosition.y() + (1 - 2*G4UniformRand())*fRasterAmplitude_vertical,
    fTargetPosition.z() 
  ); 

  //G4cout << "Updating generator..." << G4endl; 
  fEnergyAngleGenerator->Update(); 

  fParticleGun->SetParticleEnergy( fEnergyAngleGenerator->GetEnergy() );

  double cos = fEnergyAngleGenerator->GetCosTheta(); 

  G4ThreeVector momentum_direction(
    std::sqrt( 1. - cos*cos ), 
    0., 
    cos 
  );

  //rotate phi by a random angle 
  momentum_direction.rotateZ( (G4UniformRand()*2. - 1.)*CLHEP::pi/2. );

  //flip the direction if it's the RHRS 
  if (RunParameters::Instance()->Is_RHRS()) momentum_direction[1] *= -1.; 

  fParticleGun->SetParticleMomentumDirection(momentum_direction);  
  
  fParticleGun->GeneratePrimaryVertex(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
