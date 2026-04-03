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
/// \file DetectorConstruction.hh
/// \brief Definition of the B1::DetectorConstruction class

#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "UserMessenger.hh"
#include "G4MultiUnion.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

namespace B1
{

  constexpr G4double inch = 2.54*cm;

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

    void Set_BuildSieve(G4bool _val) { f_buildSieve=_val; }

    /// @return full width of sieve
    static G4double Sieve_dx() { return 4.250 * inch; }

    /// @return full height of sieve
    static G4double Sieve_dy() { return 3.250 * inch; }

    /// @return thickness of sieve
    static G4double Sieve_dz() { return 0.500 * inch; }

  protected:
    G4LogicalVolume* fScoringVolume = nullptr;
    
    /// @brief Generate multi-union solid of all sieve holes
    /// @param is_RHRS true = RHRS, false = LHRS
    /// @return A multi-union solid of all sieve-holes for the arm
    G4MultiUnion* Generate_sieveHoles_solid(const bool is_RHRS); 

    G4bool f_is_RHRS{false}; 

    /// True if the sieve-face should be constructed, false otherwise 
    G4bool f_buildSieve{true}; 

    //DetectorMessenger* fMessenger; 
    //this class' messenger 
    UserMessenger<DetectorConstruction> *fMessenger; 
};

}  // namespace B1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
