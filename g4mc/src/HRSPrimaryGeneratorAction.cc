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
/// \file HRSPrimaryGeneratorAction.cc
/// \brief Implementation of the B1::HRSPrimaryGeneratorAction class

#include "HRSPrimaryGeneratorAction.hh"
#include "RunParameters.hh"
#include "ApexTargetGeometry.hh"
#include "TargetMode.hh"
#include "EventWriter.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "TString.h"

#include <cmath> 
#include <array> 

namespace
{ 
  // maximum value of phi (azimuth around beam) from the horizontal of a generated lepton 
  constexpr G4double kPhi_max = 90.*deg; 

  // minimum angle between the beam and a generated lepton 
  constexpr G4double kTheta_min = 0.045*rad; 

  // x-position of vertical-wire targets (mm)
  constexpr G4double Vwire_X[] = {
    -3.23*mm,
    -0.72*mm,
    +1.73*mm
  };    

  // z-position of vertical-wire targets (mm)
  constexpr G4double Vwire_Z[] = {
    -196.20*mm,
      +3.80*mm,
    +203.80*mm
  };    

  //electron mass 
  constexpr double m_e = 0.501*MeV; 

  inline double RandRange(double x0, double x1) { return x0 + (x1 - x0)*G4UniformRand(); }

} // namespace


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HRSPrimaryGeneratorAction::HRSPrimaryGeneratorAction()
{
  fMessenger = new UserMessenger<HRSPrimaryGeneratorAction>(this); 

  G4String command_prefix = "/generator/"; 

  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  //get some information about this run 
  auto& run_params = RunParameters::Instance(); 

  f_sieve_x_min = run_params.Get_SieveXLow(); 
  f_sieve_x_max = run_params.Get_SieveXHigh(); 
  f_sieve_y_min = run_params.Get_SieveYLow(); 
  f_sieve_y_max = run_params.Get_SieveYHigh(); 

  f_momentum_min = run_params.Get_MomentumLow(); 
  f_momentum_max = run_params.Get_MomentumHigh();

  fTargetMode = run_params.Get_TargetMode(); 
  
  fArmMode = run_params.Get_ArmMode(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HRSPrimaryGeneratorAction::~HRSPrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HRSPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  using namespace std; 
  
  const auto& run_params = RunParameters::Instance(); 
  
  //pick a react vertex
  G4ThreeVector vertex, R_momentum, L_momentum;  

  // target code. 
  // 0 - no target specified. 
  // W - production foils 
  // V - vertical wires
  // H - horizontal wires 
  // O - optics foils
  char target_code; 

  // index of target. ex;  'V1' would have code 'V' and num '1'. 
  int target_num; 

  auto& event_writer = EventWriter::Instance(); 

  //choose a vertex 

  switch (fTargetMode)
  {
    case TargetMode::kProduction : 
      //pick a random production foil 
      int i_foil; 
      vertex = G4ThreeVector(
        (1. - 2.*G4UniformRand())*f_horizontal_raster_amplitude/2.,
        (1. - 2.*G4UniformRand())*f_vertical_raster_amplitude/2.,
        RandFoilZ(i_foil)
      ); 
      target_code ='W'; 
      target_num  =i_foil+1; 
      break; 

    case TargetMode::kV1 :
      target_code ='V'; 
      target_num  =1; 
      vertex = G4ThreeVector(
        Vwire_X[0],
        (1. - 2.*G4UniformRand())*f_vertical_raster_amplitude/2.,
        Vwire_Z[0]
      ); 
      break; 
    
    case TargetMode::kV2 :
      target_code ='V'; 
      target_num  =2; 
      vertex = G4ThreeVector(
        Vwire_X[1],
        (1. - 2.*G4UniformRand())*f_vertical_raster_amplitude/2.,
        Vwire_Z[1]
      ); 
      break; 
    
    case TargetMode::kV3 :
      target_code ='V'; 
      target_num  =3; 
      vertex = G4ThreeVector(
        Vwire_X[2],
        (1. - 2.*G4UniformRand())*f_vertical_raster_amplitude/2.,
        Vwire_Z[2]
      );
      break; 
  
    default : 
      G4Exception(
        Form("HRSPrimaryGeneratorAction::%s",__func__),
        "Unsupported target mode", 
        G4ExceptionSeverity::RunMustBeAborted, 
        Form("Target mode of Primary Generator is unsupported: '%s'",TargetMode::GetName(fTargetMode))
      );
      return; 
  }

  auto& event_writer = EventWriter::Instance(); 

  //take a position on the sieve-face, and convert it to hall coordinates (displaced from the apex scat. chamber center)
  auto sieve_pos_to_HCS = [](double x_sv, double y_sv, ArmMode::EMode mode) 
  { 
    bool is_RHRS = (mode == ArmMode::kRHRS); 

    G4ThreeVector pos_on_sieve = 
      G4ThreeVector( x_sv, y_sv, 0. ) + 
      ApexTargetGeometry::Get_sieve_pos(is_RHRS); 

      pos_on_sieve.rotateZ( -CLHEP::pi/2. );
      pos_on_sieve.rotateY( ApexTargetGeometry::Get_sieve_angle(is_RHRS) );

      return pos_on_sieve; 
  }; 
  
  //set the particle gun's position
  fParticleGun->SetParticlePosition(vertex); 

  // right-arm (positron)
  if (RHRS_is_active()) {
    
    //pick a random spot on the face of the sieve. 
    double x_sv = RandRange( f_sieve_x_min, f_sieve_x_max );
    double y_sv = RandRange( f_sieve_y_min, f_sieve_y_max );  
    
    //now, rotate this intercept back to sieve-coordinates
    auto sieve_intercept = sieve_pos_to_HCS(x_sv, y_sv, ArmMode::kRHRS); 
    
    double momentum_mag = RandRange( f_momentum_min, f_momentum_max );

    R_momentum = (sieve_intercept - vertex).unit() * momentum_mag; 

    //update data for the track data struct 
    auto& track_data = event_writer.GetTrackData(ArmMode::kRHRS_bool); 

    track_data.status         = TrackData_t::kAlive; 
    track_data.target_code    = target_code; 
    track_data.target_num     = target_num; 
    track_data.particle_type  = TrackData_t::kPositron;

    fParticleGun->SetParticleEnergy(std::sqrt( R_momentum.mag2() + m_e*m_e )); 
    fParticleGun->SetParticleMomentumDirection(R_momentum); 
    fParticleGun->SetParticleDefinition(G4Positron::Positron()); 
    fParticleGun->GeneratePrimaryVertex(event); 
  }

  // right-arm (positron) ----------------------------------------------------------------------
  if (LHRS_is_active()) {
    
    //pick a random spot on the face of the sieve. 
    double x_sv = RandRange( f_sieve_x_min, f_sieve_x_max );
    //here, because of the mirror symmetry between both arms, we have to swap min <=> max. 
    double y_sv = RandRange( -f_sieve_y_max, -f_sieve_y_min ); 
    
    //now, rotate this intercept back to sieve-coordinates
    auto sieve_intercept = sieve_pos_to_HCS(x_sv, y_sv, ArmMode::kLHRS); 
    
    double momentum_mag = RandRange( f_momentum_min, f_momentum_max );

    L_momentum = (sieve_intercept - vertex).unit() * momentum_mag; 

    //update data for the track data struct 
    auto& track_data = event_writer.GetTrackData(ArmMode::kLHRS_bool); 

    track_data.status         = TrackData_t::kAlive; 
    track_data.target_code    = target_code; 
    track_data.target_num     = target_num; 
    track_data.particle_type  = TrackData_t::kElectron;

    fParticleGun->SetParticleEnergy(std::sqrt( L_momentum.mag2() + m_e*m_e )); 
    fParticleGun->SetParticleMomentumDirection(L_momentum); 
    fParticleGun->SetParticleDefinition(G4Electron::Electron()); 
    fParticleGun->GeneratePrimaryVertex(event); 
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
