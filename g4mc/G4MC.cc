// ********************************************************************
//
// $Id: HRSMC.cc,v 1.05, 2013/10/06  HRS Exp $
// --------------------------------------------------------------

#include <stdlib.h>
#include <iostream>
//#include "HRSRootTree.hh"
#include "G4RunManager.hh"
#include "QGSP_BERT.hh"

#include "HRSUIExecutive.hh"
#include "G4UImanager.hh"

#include "HRSPhysicsList.hh"		//Jixie's physics model
#include "PhysicsList.hh"

#include "HRSMaterial.hh"
#include "HRSDetectorConstruction.hh"
#include "HRSPrimaryGeneratorAction.hh"

#include "HRSRunAction.hh"
#include "HRSEventAction.hh"
#include "HRSTrackingAction.hh"
#include "HRSSteppingAction.hh"
#include "HRSSteppingVerbose.hh"

#include "UsageManager.hh"
#include "HRSGlobal.hh"

#include "TROOT.h"

//#include "HRSRand.hh"
#include "TTimeStamp.h"   //from root, can provide ns
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include "TFileHandler.hh"

#include <memory>

#include <dlfcn.h>

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "G4String.hh"
#include <string>

using namespace std; 

//global variables
UsageManager* gConfig; //will be initialized in main() before all other classes
//HRSRootTree* gHRSTree=0; //will be initialized in main() after G4RunManager start up before user-action start

TFileHandler* fOutFile; 


int main(int argc,char** argv)
{
  const char *const here = "G4MC::main"; 
  
  
  //  dlopen("libpython2.7.so", RTLD_NOW | RTLD_NOLOAD | RTLD_GLOBAL);
  ////////////////////////////////////////////////////////////////////
  //take care of the arguments
  ////////////////////////////////////////////////////////////////////
  string LogFileName("G2P_G4Sim_Log.txt");
  gConfig= new UsageManager("HRSUsage.ini",argc,argv,LogFileName.c_str());
  gConfig->ReadFile("Detector.ini");
  //gConfig->PrintOpt();
  //gConfig->PrintParamMap(); 
  //  if (seed==0)

  //root file to use
  //gConfig
  
  unsigned Seed;
  TTimeStamp pTS;
  Seed=unsigned(time(0))+pTS.GetNanoSec();
  cout<<"seed = "<<Seed<<endl;
  srand(Seed);
  cout<<" seed rand = "<<rand()%1000000<<endl;
  G4long myseed = Seed;
  CLHEP::HepRandom::setTheSeed(myseed);
  
  //  else Seed=seed;
  
  //printf("\n***HRSRand::HRSRand(): set the random seed to %d ***\n",Seed);
  
  //  mSeedLR=rand();
  
  ////////////////////////////////////////////////////////////////////
  //my Stepping verbose output class
  ////////////////////////////////////////////////////////////////////
  G4VSteppingVerbose::SetInstance(new HRSSteppingVerbose);
  
  ////////////////////////////////////////////////////////////////////
  // RunManager construction
  ////////////////////////////////////////////////////////////////////
  auto runManager = new G4RunManager;
  
  
  ////////////////////////////////////////////////////////////////////
  // Visualization manager construction	
  ////////////////////////////////////////////////////////////////////
#ifdef G4VIS_USE
  auto visManager = new G4VisExecutive;
  //default verbose level= warning, 
  visManager->Initialize();
  
  cout << "HERE ~~~~~~~~~~~~~ vis manager initialized" << endl; 
#endif
  
  ////////////////////////////////////////////////////////////////////
  // mandatory user initialization classes
  ////////////////////////////////////////////////////////////////////
  //initial materials
  auto pMaterialManager = new HRSMaterial;
  
  runManager->SetUserInitialization(new HRSDetectorConstruction);
  
  //inital physics list
  int pUseJixieModel=0;
  
  gConfig->GetArgument("UseJixieModel",pUseJixieModel);
  
  if(pUseJixieModel==1) {
    cout<<"using Jixie physics list"<<endl;
    //Use Jixie's model
    runManager->SetUserInitialization(new HRSPhysicsList);
  }
  
  /*
    else
gConfig->GetArgument("NoSecondary",iNoSecondary);    {
    //or use this one which is from example/extended/hadon/had01
    G4String pPhysicsModel=gConfig->GetArgument("PhysicsModel");
    runManager->SetUserInitialization(new PhysicsList(pPhysicsModel));
    }
  */

  if(pUseJixieModel==-1) {
    cout<<"using default QGSP_BERT"<<endl;
    G4VUserPhysicsList* physics_l = new QGSP_BERT();
    runManager->SetUserInitialization(physics_l);
  }

  // initialize Geant4 kernel
  runManager->Initialize();
  // mandatory user action class
  auto hrs_primary_generator_action = new HRSPrimaryGeneratorAction; 
  
  runManager->SetUserAction(hrs_primary_generator_action);
  
  
  // The analysis manager
  int pRunNumber=1;
  gConfig->GetArgument("RunNumber",pRunNumber);
  /*
    gHRSTree=new HRSRootTree(pRunNumber); 
  */
  //if(pRunNumber<0) gHRSTree->SetRunNumber(abs(pRunNumber)); 

  ////////////////////////////////////////////////////////////////////
  // optional user action classes, 
  ////////////////////////////////////////////////////////////////////
  // runManager will delete them during closing
  runManager->SetUserAction(new HRSRunAction);
  runManager->SetUserAction(new HRSEventAction);
  runManager->SetUserAction(new HRSTrackingAction);
  //auto stepping_action = new HRSSteppingAction;
  runManager->SetUserAction(new HRSSteppingAction);
  //runManager->SetUserAction(stepping_action);

  ////////////////////////////////////////////////////////////////////
  // User interactions
  //////////////////////////////////////////////////////////////////// 
  // execute macro files if provided through arguments
  // Define (G)UI for interactive mode
  // The UImanager will be handled by Ranmanager,no need to delete it manually 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  char cmd[255],tmpstr[255];
  //take care of cmd files then the mcin file
  int pNofMac=2;
  //gConfig->GetArgument("NofMac",pNofMac); 

  

  /*std::string MacFile;
  for(int i=0;i<pNofMac;i++) {
    sprintf(tmpstr,"MacFile%d",i+1);
    MacFile=gConfig->GetArgument(tmpstr);
    if(CheckFile(MacFile)) {
      sprintf(cmd,"/control/execute %s",MacFile.c_str()); 
      cout << "Trying to execute cmd: ("<<i<<") " << cmd << endl;
      UImanager->ApplyCommand(cmd);
    }
    }*/

  //Hard-coding this bit. the arguments are going to just be first two macros fed
  // to the executable. we expect 2 macros.
  
  //this macro configures the detectors / file readers
  string macro_name = gConfig->GetArgument("MacFile1"); 
  
  printf("<%s>: Attempting to execute macro '%s'...\n", here, macro_name.data());   
  
  sprintf(cmd,"/control/execute %s",macro_name.data()); 

  UImanager->ApplyCommand(cmd);

  
  //auto out_file = unique_ptr<TFileHandler>(new TFileHandler("out_Q1.root"));
  string outfile_path = hrs_primary_generator_action->Get_outfile_path().data(); 
  
  fOutFile = new TFileHandler( outfile_path.data(),
			       hrs_primary_generator_action->Is_RHRS() ); 
  
  printf("<%s>: Outfile path is: %s\n", here, outfile_path.data()); 
  
  printf("Gun Z Lo/Hi :[% 5.f, % 5.f]\n",
	 hrs_primary_generator_action->GetGunZLow(),
	 hrs_primary_generator_action->GetGunZHigh()); 

  
  macro_name = gConfig->GetArgument("MacFile2"); 
  
  printf("<%s>: Attempting to execute macro '%s'...\n", here, macro_name.data());   
  
  sprintf(cmd,"/control/execute %s",macro_name.data()); 
  
  UImanager->ApplyCommand(cmd);

  
  
  
  //By Jixie: Might change the trigger for this part later
  //gConfig->GetArgument("UseRootNtuple",pUseRootNtuple); 
  //if(UseRootNtuple) //root ntuple mode	
  string pPrimaryEngine1=gConfig->GetArgument("PrimaryEngine1"); 
  if(pPrimaryEngine1=="RootNtuple") //root ntuple mode	
    {		
      int pTrigNum=-1;
      gConfig->GetArgument("TrigNum1",pTrigNum);
      if(pTrigNum<=0) pTrigNum=9999999;
      sprintf(cmd,"/run/beamOn %d",pTrigNum); 
      G4cout<<"Trying to execute cmd: "<<cmd<<G4endl;
      UImanager->ApplyCommand(cmd);
    }
  ////////////////////////////////////////////////////////////////////
  // interaction mode
  ////////////////////////////////////////////////////////////////////

  
  

  
  int pInteractiveMode=0;
  gConfig->GetArgument("InteractiveMode",pInteractiveMode); 
  if(pInteractiveMode)
    {
      int pUseGui=0;
      gConfig->GetArgument("UseGui",pUseGui);
      HRSUIExecutive *pUI = new HRSUIExecutive(argc, argv, pUseGui);
      if(pUI->IsGUI())
	{
	  //put your extra cmd here;
	  if(CheckFile("gui.mac"))
	    {
	      UImanager->ApplyCommand("/control/execute gui.mac");
	    }
	}
      pUI->SessionStart();
    }
  ////////////////////////////////////////////////////////////////////
  //free the memory
  ////////////////////////////////////////////////////////////////////
  //	if(gHRSTree) delete gHRSTree;

#ifdef G4VIS_USE
  delete visManager;
  //G4cout << "Vis manager deleted..." << G4endl;
#endif

  // Free the memory: user actions, physics_list and detector_description are
  //                  owned and deleted by the run manager, so they should not
  //                  be deleted in the main() program !

  delete pMaterialManager;

  delete runManager;
  //G4cout << "Run manager deleted..." << G4endl;

  fOutFile->CloseFile(); 
  delete fOutFile; 
  
  return 0;
}

