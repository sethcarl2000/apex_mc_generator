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
/// \file SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "RunParameters.hh"
#include "EventAction.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4ThreeVector.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Step.hh"
#include "G4Track.hh"

#include <math.h>

namespace {
  //the minimum angle between the beam and an electron/positron that may make it into the detector
  const G4double max_cosine = std::cos( 0.045 ); 
}

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction) : fEventAction(eventAction) {

  // tan(theta_min)
  fMinTanTheta = std::tan( RunParameters::Instance()->GetMinAngle() ); 

  fMomentum_min = RunParameters::Instance()->GetMomentum_min(); 
  fMomentum_max = RunParameters::Instance()->GetMomentum_max(); 

  fz_target_back = RunParameters::Instance()->GetTargetThickness()/2; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......s

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  //ask the detector construction class for a ptr to the scoring volume
  if (!fScoringVolume) {
    const auto detConstruction = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();

    //ask the detector construction which arm we're using 
    f_is_RHRS = RunParameters::Instance()->Is_RHRS(); 
  }

  // get volume of the current step
  G4LogicalVolume* volume =
    step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  // check if we are in scoring volume
  if (volume != fScoringVolume) return;

  G4Track* track = step->GetTrack(); 

  // kills current track 
  auto kill_track = [track]() { track->SetTrackStatus(fStopAndKill); };
  
  G4int charge{0}; 
  auto particle_def = track->GetParticleDefinition(); 

  //reject any particles which aren't positrons or electrons
  if (particle_def == G4Electron::Electron()) {
    charge = -1; 
  } else {
    if (particle_def == G4Positron::Positron()) {
      charge = +1;
    } else { kill_track(); return; }
  }
  
  const G4ThreeVector& p_track = track->GetMomentum();  

  //check if the track is inside the momentum acceptance
  const RunParameters* run_params = RunParameters::Instance(); 

  G4double p_track_mag = p_track.mag(); 
  
  //only save leptons inside our momentum acceptance 
  if (p_track_mag < fMomentum_min || p_track_mag > fMomentum_max ) { kill_track(); return; }  

  //only save leptons whose angle with the beam is beyond a certain value
  G4double tan_theta = sqrt( p_track[0]*p_track[0] + p_track[1]*p_track[1] ) / p_track[2]; 

  if (p_track[2] < 0 || tan_theta < fMinTanTheta) { kill_track(); return; }
      
  //get the analysis manager, and fill out relevant information. 
  auto analysisManager = G4AnalysisManager::Instance(); 

  int i_col=0; 

  //save the momentum 
  analysisManager->FillNtupleDColumn(i_col++, p_track.x());
  analysisManager->FillNtupleDColumn(i_col++, p_track.y());
  analysisManager->FillNtupleDColumn(i_col++, p_track.z());

  analysisManager->FillNtupleDColumn(i_col++, p_track_mag); //momentum magnitude
  analysisManager->FillNtupleDColumn(i_col++, std::atan(tan_theta)); //angle between track and z-axis
  
  const G4ThreeVector& r_track = track->GetPosition(); 

  //save the position. project the particle onto the back of the target 
  double dxdz = p_track.x() / p_track.z(); 
  double dydz = p_track.y() / p_track.z(); 
  
  analysisManager->FillNtupleDColumn(i_col++, r_track.x() - dxdz*(r_track.z() - fz_target_back));
  analysisManager->FillNtupleDColumn(i_col++, r_track.y() - dydz*(r_track.z() - fz_target_back));
  analysisManager->FillNtupleDColumn(i_col++, fz_target_back);

  //save the charge to identify what kind of particle it is 
  analysisManager->FillNtupleIColumn(i_col++, charge);

  //save this as a distinct event 
  analysisManager->AddNtupleRow(); 

  //now, get rid of this track, to avoid saving it twice. 
  kill_track(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
