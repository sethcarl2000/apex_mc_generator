// ********************************************************************
//
// $Id: HRSCalorimeterSD.cc,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//
#include "HRSCalorimeterSD.hh"
#include "HRSCalorimeterHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//By Jixie: The calorimeter SD is only goof for EC. 
//It will not pay attention to the time window. 
// For a calorimeter, I assume each sensitive detector is just one block
// It will create only one calorimeter hit for each SD in each event 
// (one hit one event)
// For this calorimeter hit, it will accumulate all deposited energy from all 
// sources in this event. It could include energy from multiple tracks.
// Because of this design, I'd better to shoot only one primary particle
// to EC in each event 
// This hit will record the inpos and inmom from the initial interatction 
// by one track and  will record the time, outpos and outmom from the last 
// interaction os the same track (usually the first track firing it)
// In the event action, these hits will be stored into the root ntuple


HRSCalorimeterSD::HRSCalorimeterSD(G4String name):G4VSensitiveDetector(name)
{
	G4String HCname;
	collectionName.insert(HCname="CalorimeterColl");
	HCID = -1;
	//The given 'name' will be stored into G4VSensitiveDetector::SensitiveDetectorName
	//you can give '/RTPCGEM' or 'RTPCGEM', the Geant4 system will both set 'RTPCGEM'
	//as the SensitiveDetectorName and '/RTPCGEM' as fullPathName
	//
	//For current hit colection, the collectionName is "SensitiveDetectorName/HCname"
	//To access this HC in event action, one need to get HCID through
	//G4HCtable::GetCollectionID(G4String collectionName) 

	hitsCollection = 0;
	verboseLevel = 0; //defined in G4VSensitiveDetector.hh
}

HRSCalorimeterSD::~HRSCalorimeterSD()
{
	;
}

void HRSCalorimeterSD::Initialize(G4HCofThisEvent* HCE)
{
	hitsCollection = new HRSCalorimeterHitsCollection(SensitiveDetectorName,collectionName[0]);
	if(HCID<0)
	{ 
		HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); 
	}
	HCE->AddHitsCollection(HCID,hitsCollection);
}


G4bool HRSCalorimeterSD::ProcessHits(G4Step* aStep,G4TouchableHistory* /*aROHist*/)
{	
	if(!hitsCollection) return false;

	G4double edep = aStep->GetTotalEnergyDeposit();
	if(edep/keV<=0.)  return true;


	G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
	G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
	G4int copyNo = theTouchable->GetVolume()->GetCopyNo();
	G4int trackId = aStep->GetTrack()->GetTrackID();
	G4double hitTime = preStepPoint->GetGlobalTime();
	
	G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
	G4ThreeVector outpos = postStepPoint->GetPosition();       
	G4ThreeVector outmom = postStepPoint->GetMomentum();
	G4String physName = theTouchable->GetVolume()->GetName();

	//By Jixie; it turns out that this is always zero due to the fact that G4 physics process never set 
	//a value to it
	G4double edep_NonIon = aStep->GetNonIonizingEnergyDeposit();
	G4String process=postStepPoint->GetProcessDefinedStep()->GetProcessName();
	if(process=="Transportation" || process=="msc") edep_NonIon=edep;	
	//if(edep_NonIon) G4cout<<"edep_NonIon="<<edep_NonIon<<G4endl;


	// check if this finger already has a hit
	// generally a hit is the total energy deposit in a sensitive detector within a time window
	// it does not matter which particle deposite the energy
	// For a calorimeter, I assume each sensitive detector is just one block
	// Therefore I will create only one calorimeter hit for each SD
	// For this calorimeter hit, I will accumulate all deposited energy from all 
	// sources in this event. It could include energy from multiple tracks.
	// Because of this design, I'd better to shoot only one primary particle
	// to EC in each event 

	HRSCalorimeterHit* aHit = 0;
	for(G4int i=hitsCollection->entries()-1;i>=0;i--)
	{		
		if( (*hitsCollection)[i]->GetPhysV()->GetName()==physName && 
			(*hitsCollection)[i]->GetId()==copyNo)
		{
			//found an exist hit 
			if(verboseLevel >= 3)
			{
				G4cout<<"found an exist Calorimeter hit: PhysVol = "<<physName
					<<" copyNo = "<<theTouchable->GetVolume()->GetCopyNo()<<G4endl;
			}
			aHit = (*hitsCollection)[i];
			break;
		}
	}
	

	// if it has, then take the earlier time, accumulated the deposited energy
	// if not, create a new hit and set it to the collection
	if(aHit)
	{
		if(verboseLevel >= 4)
		{
			G4cout<<"<<<add energy into an exist hit: position="<<outpos<<"  edep="<<edep<<G4endl;
		}
		aHit->AddEdep(edep); 
		aHit->AddNonIonEdep(edep_NonIon);
		if(aHit->GetTrackId()==trackId)
		{
			if(aHit->GetTime() > hitTime) aHit->SetTime(hitTime); 
			aHit->SetOutPos(outpos);
			aHit->SetOutMom(outmom);
		}
	}
	else
	{
		//this is a new hit
		aHit = new HRSCalorimeterHit(copyNo,hitTime);

		G4int parentTrackId=aStep->GetTrack()->GetParentID();
		G4ThreeVector inpos = preStepPoint->GetPosition();       
		G4ThreeVector inmom = preStepPoint->GetMomentum();

		G4int pdgid=aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
		if(pdgid>3000)
		{
			//G4cout<<"Wrong PDGCCoding pdgid="<<pdgid
			//	<<" name="<<aStep->GetTrack()->GetParticleDefinition()->GetParticleName()<<G4endl;
			pdgid=(pdgid%1000000)/10;
		}

		aHit->SetPhysV(theTouchable->GetVolume());
		aHit->SetEdep(edep);
		aHit->SetNonIonEdep(edep_NonIon);
		aHit->SetTrackId(trackId);
		aHit->SetOutPos(outpos);
		aHit->SetOutMom(outmom);
		aHit->SetInPos(inpos);
		aHit->SetInMom(inmom);
		aHit->SetParentTrackId(parentTrackId);
		aHit->SetPdgid(pdgid);

		hitsCollection->insert(aHit);

		if(verboseLevel >= 2)
		{
			G4cout<<"<<<Create a new Hit of type <" << SensitiveDetectorName << "/"<<collectionName[0]
			<<" in "<<theTouchable->GetVolume()->GetName()<<"["<<copyNo<<"] \n"
				<<" at "<<inpos<<", by track "<< trackId <<", momentum="<<inmom<<G4endl;
		}
	}

	return true;
}

void HRSCalorimeterSD::EndOfEvent(G4HCofThisEvent* HCE)
{
	if(!HCE) return;  //no hits collection found 
	//the above line just to avoid warning of not use HCE

	int nhitC = hitsCollection->GetSize();
	if(!nhitC) return;

	if(verboseLevel >= 2 && nhitC)
	{
		G4cout<<"<<<End of Hit Collections <" << SensitiveDetectorName << "/"<<collectionName[0]
		<<">: " << nhitC << " hits." << G4endl; 
		for (int i=0; i<nhitC; i++)
		{
			(*hitsCollection)[i]->Print();
		}
	}

	return;
}

