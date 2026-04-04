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
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "ApexTargetGeometry.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/RotationZ.h"
#include "CLHEP/Vector/RotationX.h"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4MultiUnion.hh"
#include "G4Trd.hh"
#include "RunParameters.hh"

#include <stdio.h> 

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction()
{
  //fMessenger = new DetectorMessenger(this); 
  fMessenger = new UserMessenger<DetectorConstruction>(this); 
  
  G4String command_prefix = "/detector/"; 

  fMessenger->AddCommand_bool(
    command_prefix + "build_sieve",
    "build_sieve",
    &DetectorConstruction::Set_BuildSieve,
    true,
    "if 'true', then a solid tungsten sieve will be constructed. if 'false', no solid sieve will be constructed"
  ); 
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{
  delete fMessenger; 
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



G4VPhysicalVolume* DetectorConstruction::Construct()
{
  const RunParameters* run_params = RunParameters::Instance(); 

  const bool is_RHRS = run_params->Is_RHRS(); 
  
  // Get nist material manager
  // 
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  // 
  G4double env_sizeXY = 10 * cm, env_sizeZ = 10 * cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_Galactic"); /// this represents a high-vaccuum environment 

  // Option to switch on/off checking of volumes overlaps
  // 
  G4bool checkOverlaps = true;
  
  // 
  // Tungsten sieve face 
  // 
  G4Material* tungsten_mat = nist->FindOrBuildMaterial("G4_W"); //tungsten material
  
  //
  // World
  //
  G4double world_sizeXY =  50*cm;
  G4double world_sizeZ  = 200*cm;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

  auto solidWorld = new G4Box("World",  // its name:
    // its dimensions
    world_sizeXY/2., world_sizeXY/2., world_sizeZ/2.
  ); 

  auto logic_World = new G4LogicalVolume(
    solidWorld,  // its solid
    world_mat,  // its material
    "World");  // its name

  auto phys_World = new G4PVPlacement(
    nullptr,  // no rotation
    G4ThreeVector(),  // at (0,0,0)
    logic_World,  // its logical volume
    "World",  // its name
    nullptr,  // its mother  volume
    false,  // no boolean operation
    0,  // copy number
    checkOverlaps  // overlaps checking
  ); 
  
  /// dimensions of the sieve 
  const G4double sieve_dx = Sieve_dx();
  const G4double sieve_dy = Sieve_dy(); 
  const G4double sieve_dz = Sieve_dz(); 

  const G4double sieve_y_angle = ApexTargetGeometry::Get_sieve_angle(is_RHRS); 

  /// Scoring volume -------------------------------------
  //
  // This volume will monitor the tracks, and record ones we want to keep (and kill the ones we don't)
  const G4double scoring_volume_thickness = 2.5*mm; 
  auto solid_scoringVolume = new G4Box(
    "solid_scoringVolume",
    20*cm, 
    20*cm, 
    scoring_volume_thickness/2.
  ); 
  auto logic_scoringVolume = new G4LogicalVolume(
    solid_scoringVolume,
    world_mat,
    "logic_scoringVolume"
  ); 

  //only show the scoring box if the solid sieve was not built 
  new G4PVPlacement(
    nullptr, 
    G4ThreeVector(0., 0., run_params->GetTargetThickness()/2 + 2*mm), // let's leave a 2 cm gap between the target and the sieve
    logic_scoringVolume, 
    "Scoring Volume", 
    logic_World, 
    true, 
    0, 
    checkOverlaps
  );
  fScoringVolume = logic_scoringVolume; 



  G4VSolid *solid_target = nullptr; 
  G4RotationMatrix *rotation_target = nullptr; 

  solid_target = new G4Box(
    "solid_Target", 
    10*cm/2.,
    10*cm/2.,
    run_params->GetTargetThickness()/2.
  ); 

  //rotation_target = new G4RotationMatrix( CLHEP::HepRotationX(CLHEP::pi/2.) ); 
  
  auto logic_target = new G4LogicalVolume(
    solid_target, 
    tungsten_mat, 
    G4String("logic_target")
  ); 

  new G4PVPlacement(
    nullptr, //rotation_target,
    G4ThreeVector(0, 0, 0),
    logic_target, 
    "Target", 
    logic_World, 
    false, 
    0, 
    false
  );

  //
  //  always return the physical World
  //
  return phys_World;
}

G4MultiUnion* DetectorConstruction::Generate_sieveHoles_solid(const bool is_RHRS)
{     
  const G4double sieve_dz = 0.500 * inch; 

  const G4double holeRad_small = 0.055/2. * inch; 
  const G4double holeRad_back  = 0.075/2. * inch; 
  const G4double holeRad_big   = 0.106/2. * inch; 

  /// A 'small' sieve-hole with a constant bore diameter
  auto solid_sieveHole_small = new G4Tubs(
    "solid_sieveHole_small", 
    0., holeRad_small, 
    1.2*sieve_dz,
    0., 2.*CLHEP::pi 
  );

  /// A 'big' sieve-hole with a constnat bore diameter
  auto solid_sieveHole_big = new G4Tubs(
    "solid_sieveHole_small", 
    0., holeRad_big, 
    1.2*sieve_dz,
    0., 2.*CLHEP::pi 
  );

  /// A sieve-hole with a wider diameter at the back 
  auto solid_sieveHole_wideBack = new G4UnionSolid(
    "soild_sieveHole_wideBack",
    solid_sieveHole_small, 
    new G4Tubs(               //this is the wider-radius downstream part of the hole 
      "solid_sieveHole_back",
      0., holeRad_back, 
      1.2*sieve_dz/2., 
      0., 2.*CLHEP::pi
    ),
    G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., +1.2*sieve_dz/2.))
  );

  //now, let's put them all in one place 
  auto solid_allSieveHoles = new G4MultiUnion("solid_allSieveHoles"); 

  const double dx = 0.230 * inch; //all units in meters

  //the sign is different here, because the L & R-sieves are mirror-images of
  // of one another
  const double dy = is_RHRS ? +0.190 * inch : -0.190 * inch;
  
  //the first row is 8 rows above (-x) the center hole
  const double x0 = -dx * 8.; 
  //the first column is 7-spaces +y from the center
  const double y0 =  dy * 7.; 

  const int nRows=17; 
  for (int row=0; row<nRows; row++) { 

    //17 rows in total, numbered 0-16.
    // row 0 is the highest in HCS, so lowest-X in TCS.
    // Recall that in TCS, the central 'big hole' is the origin of x & y.
    
    //even rows are shifted in +y half a col-spacing
    bool evenRow = (row % 2==0); 
    int nCols;
    if (evenRow) { nCols = 15; }
    else  {
        if (row==1 || row==15) {nCols = 12;} //each odd row has a 'gap' missing, except these two
        else                   {nCols = 11;}
    } 
    
    for (int col=0; col<nCols; col++) { 
    
      //check wheter this hole is a 'big' hole
      bool bigHole ((row==8 && col==7) || (row==12 && col==3)); 

      //get this hole's x-position
      double x = x0 + ((double)row)*dx;
      double y = y0 - ((double)col)*dy;

      if (!evenRow) { 
          y += -dy/2.; //odd-holes are shifted in y a bit
          //skip over the gap which happens in some rows, but not these. 
          if ( (row!=1 && row!=15) && col>5 ) y += -dy; 
      }

      //find out what kind of sieve hole this is 
      G4VSolid* hole_solid; 
      
      if (bigHole) {
        hole_solid = solid_sieveHole_big;
      } else {
        if ((row <= 2 || row >= 14) && col <= 12) {
          hole_solid = solid_sieveHole_wideBack; 
        } else {
          hole_solid = solid_sieveHole_small; 
        }
      } 

      solid_allSieveHoles->AddNode(
        hole_solid, 
        G4Transform3D(G4RotationMatrix(), G4ThreeVector(x, y, 0.))
      );

    }//for (int col=0; col<nCols; col++) 

  }//for (int row=0; row<nRows; row++) 

  solid_allSieveHoles->Voxelize(); 

  return solid_allSieveHoles; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
