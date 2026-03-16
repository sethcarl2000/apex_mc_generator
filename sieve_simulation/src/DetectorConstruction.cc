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

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4MultiUnion.hh"
#include "G4Trd.hh"

namespace {
  constexpr G4double inch = 2.54 * cm; 

}

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  const bool is_RHRS = false; 
  
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
  G4double world_sizeXY = 10. * cm;
  G4double world_sizeZ  = 10. * cm;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

  auto solidWorld = new G4Box(/*its name: */"World",  /*its dimensions: */world_sizeXY/2., world_sizeXY/2., world_sizeZ/2.); 

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
  
  /// Position of the sieve
  G4ThreeVector position_sieve(0., 0., 0.); 
  /// Rotation of the sieve 
  G4RotationMatrix *rotation_sieve = nullptr; 
  
  /// dimensions of the sieve 
  const G4double sieve_dx = 3.250 * inch; 
  const G4double sieve_dy = 4.250 * inch; 
  const G4double sieve_dz = 0.500 * inch; 

  /// Sieve containter volume 
  // Physical volume
  auto solid_SieveContainer = new G4Box(
    "sieve_container", 
    (sieve_dx + 1*cm)/2.,
    (sieve_dy + 1*cm)/2.,
    (sieve_dz + 1*cm)/2. 
  );
  // logical volume
  auto logic_sieveContainer = new G4LogicalVolume(
    solid_SieveContainer, 
    world_mat, 
    "logic_sieveContainer"
  );
  //make the sieve container invisible
  logic_sieveContainer->SetVisAttributes(G4VisAttributes::GetInvisible()); 

  //place the sieve container in the world
  new G4PVPlacement(
    rotation_sieve, 
    position_sieve, 
    logic_sieveContainer,
    "Sieve Container",
    logic_World,
    true, 
    0, 
    checkOverlaps
  ); 

  /// sieve face -----------------------------------------

  // Solid volume
  auto solid_sieveFace = new G4Box(
    "solid_sieveFace", 
    sieve_dx/2., 
    sieve_dy/2., 
    sieve_dz/2.
  );   
  // first, we create 3 different types of sieve-holes:

  auto solid_allSieveHoles = Generate_sieveHoles_solid(is_RHRS); 

  // make a 'subtraction solid' of the sieve face (sieve face with holes drilled into it )
  auto solid_sieveWithHoles = new G4SubtractionSolid(
    "solid_sieveWithHoles", 
    solid_sieveFace, 
    solid_allSieveHoles, 
    G4Transform3D(G4RotationMatrix(), G4ThreeVector())
  ); 
  //solid_sieveWithHoles->CreatePolyhedron(); 

  // logical volume 
  auto logic_sieveFace = new G4LogicalVolume(
    solid_sieveWithHoles, 
    tungsten_mat, 
    "logic_sieveFace"
  );
  logic_sieveFace->SetVisAttributes(G4VisAttributes(true)); 
  //add the fsieve face as a placement 
  new G4PVPlacement(
    nullptr,         // no rotation
    G4ThreeVector(0., 0., 0.), // placement at origin
    logic_sieveFace,
    "Sieve Face",
    logic_sieveContainer, 
    true, 
    0, 
    checkOverlaps
  ); 
  // -----------------------------------------------------

/*  
  //
  // Shape 1
  //
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4ThreeVector pos1 = G4ThreeVector(0, 2 * cm, -7 * cm);

  // Conical section shape
  G4double shape1_rmina = 0. * cm, shape1_rmaxa = 2. * cm;
  G4double shape1_rminb = 0. * cm, shape1_rmaxb = 4. * cm;
  G4double shape1_hz = 3. * cm;
  G4double shape1_phimin = 0. * deg, shape1_phimax = 360. * deg;
  auto solidShape1 = new G4Cons("Shape1", shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb,
                                shape1_hz, shape1_phimin, shape1_phimax);

  auto logicShape1 = new G4LogicalVolume(solidShape1,  // its solid
                                         shape1_mat,  // its material
                                         "Shape1");  // its name

  new G4PVPlacement(nullptr,  // no rotation
                    pos1,  // at position
                    logicShape1,  // its logical volume
                    "Shape1",  // its name
                    logicEnv,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking

  //
  // Shape 2
  //
  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  G4ThreeVector pos2 = G4ThreeVector(0, -1 * cm, 7 * cm);

  // Trapezoid shape
  G4double shape2_dxa = 12 * cm, shape2_dxb = 12 * cm;
  G4double shape2_dya = 10 * cm, shape2_dyb = 16 * cm;
  G4double shape2_dz = 6 * cm;
  auto solidShape2 =
    new G4Trd("Shape2",  // its name
              0.5 * shape2_dxa, 0.5 * shape2_dxb, 0.5 * shape2_dya, 0.5 * shape2_dyb,
              0.5 * shape2_dz);  // its size

  auto logicShape2 = new G4LogicalVolume(solidShape2,  // its solid
                                         shape2_mat,  // its material
                                         "Shape2");  // its name

  new G4PVPlacement(nullptr,  // no rotation
                    pos2,  // at position
                    logicShape2,  // its logical volume
                    "Shape2",  // its name
                    logicEnv,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking

  // Set Shape2 as scoring volume
  //
*/ 
  fScoringVolume = logic_sieveContainer;
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
