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

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction) : fEventAction(eventAction) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);

  // post step point: 
  //const G4StepPoint* postStepPoint = step->GetPostStepPoint(); 

  G4Track* track = step->GetTrack(); 

  const G4ThreeVector& p_track = track->GetMomentum();  

  // kills current track 
  auto kill_track = [track]() { track->SetTrackStatus(fStopAndKill); };

  //check if the track is inside the momentum acceptance
  const RunParameters* run_params = RunParameters::Instance(); 

  G4double p_track_mag = p_track.mag(); 
  
  //only save leptons inside our momentum acceptance 
  /*if (p_track_mag < run_params->GetMomentum_min() || 
      p_track_mag > run_params->GetMomentum_max() ) { kill_track(); return; }*/ 
      
  //only save positrons for the right arm, and electrons for the left
  /*if (f_is_RHRS) {
    if (track->GetParticleDefinition() != G4Positron::Positron() ) { kill_track(); return; }
  } else {
    if (track->GetParticleDefinition() != G4Electron::Electron() ) { kill_track(); return; }  
  }*/ 

      
  //get the analysis manager, and fill out relevant information. 
  auto analysisManager = G4AnalysisManager::Instance(); 

  int i_col=0; 

  //save the momentum 
  analysisManager->FillNtupleDColumn(i_col++, p_track.x());
  analysisManager->FillNtupleDColumn(i_col++, p_track.y());
  analysisManager->FillNtupleDColumn(i_col++, p_track.z());
  
  const G4ThreeVector& r_track = track->GetPosition(); 

  //save the position 
  analysisManager->FillNtupleDColumn(i_col++, r_track.x());
  analysisManager->FillNtupleDColumn(i_col++, r_track.y());
  analysisManager->FillNtupleDColumn(i_col++, r_track.z());

  //save this as a distinct event 
  analysisManager->AddNtupleRow(); 

  //now, get rid of this track, to avoid saving it twice. 
  kill_track(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
