// ********************************************************************
//
// $Id: HRSEMField.hh,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//

#ifndef HRSEMField_H
#define HRSEMField_H

#include "G4ThreeVector.hh"
#include "G4ElectroMagneticField.hh"
#include "BField_Helm.hh"   //Target Field class
#include "BField_Septum.hh" //Septum Field class
#include "BField_Septum_New.hh" //Septum Field class

class HRSEMFieldMessenger;
class HRSEMField : public G4ElectroMagneticField
{
public:

	HRSEMField() ;                
	~HRSEMField() ;  

	inline void GetFieldValue(const G4double Point[4], G4double *Bfield ) const;
	//  Point[4] x,y,z,time
	//  Return as Bfield[0], [1], [2] the magnetic field x, y & z components
	//   and   as Bfield[3], [4], [5] the electric field x, y & z components

	G4bool DoesFieldChangeEnergy() const { return true; }
	//  For field with an electric component this should be true
	//  For pure magnetic field this should be false
	//  Alternative: default safe implementation { return true; }

	inline void SetErDC(G4double val) { ErDC = val; }
	inline G4double GetErDC() const { return ErDC; }

	inline void SetErInner(G4double val) { ErInner = val; }
	inline G4double GetErInner() const { return ErInner; }

	inline void SetEField3V(G4ThreeVector v) { EField3V = v; bUseUniformEField=true;}
	inline G4ThreeVector GetEField3V() const { return EField3V; }

	inline void SetBField3V(G4ThreeVector v) { BField3V = v; bUseUniformBField=true;}
	inline G4ThreeVector GetBField3V() const { return BField3V; }

private:
	HRSEMFieldMessenger* messenger;

	bool bUseUniformEField;
	bool bUseUniformBField;
        G4double mLHRSMomentum;
        G4double mRHRSMomentum;
        G4double mLHRSAngle;
        G4double mRHRSAngle;
        G4int mFringeField;
        G4int Septum_On;
        G4int SeptumNew;
        G4double SeptumFieldScale;
        G4int nQ1Sos;
        BField_Septum* b_septum;
        BField_Septum_New* b_septum_new;

	BField_Helm*    mBField_Helm;
	BField_Septum*  mBField_Septum;

	G4double ErDC;
	G4double ErInner;
	G4ThreeVector EField3V;
	G4ThreeVector BField3V;
        G4double KAPPA1;
        G4double KAPPA2;
        G4double KAPPA3;
        G4double DipField;

};

#endif
