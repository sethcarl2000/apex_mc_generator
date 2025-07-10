// ********************************************************************
//
// $Id: HRSPrimaryGeneratorMessenger.cc,v 1.0, 2010/12/26  HRS Exp $
// --------------------------------------------------------------
//
#include "HRSPrimaryGeneratorMessenger.hh"
#include "HRSPrimaryGeneratorAction.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4ThreeVector.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"
#include "G4ios.hh"

using g4uic_double_unit = G4UIcmdWithADoubleAndUnit; 
using g4uic_int         = G4UIcmdWithAnInteger;


HRSPrimaryGeneratorMessenger::HRSPrimaryGeneratorMessenger(HRSPrimaryGeneratorAction * mpga)
:target(mpga)
{
  mydetDir = new G4UIdirectory("/mydet/");
  mydetDir->SetGuidance("Jixie's commands for this program only");


  isRHRS_Cmd = unique_ptr<G4UIcmdWithABool>(new G4UIcmdWithABool("/mydet/use_RHRS",this)); 

  outfileCmd = unique_ptr<G4UIcmdWithAString>(new G4UIcmdWithAString("/mydet/outfile",this)); 
  char strCmd[100],strGuid[255];

  for(int i=0;i<MaxPrimaryNum;i++) {	
    sprintf(strCmd,"/mydet/particle%d/",i+1);
    sprintf(strGuid,"Particle gun control commands for primary particle %d",i+1);
    gunDir[i] = new G4UIdirectory(strCmd);
    gunDir[i]->SetGuidance(strGuid);
      
    //q1 not the best idea
    //		sprintf(strCmd,"/mydet/particle%d/q1grad",i+1);
    //		sprintf(strGuid,"Upper limit of theta for primary particle %d, triggered by theta<0",i+1);
    //		q1grad[i] = new g4uic_double_unit(strCmd,this);
    //		q1grad[i]->SetGuidance(strGuid);
    //		q1grad[i]->SetParameterName("q1grad",false);
    //		q1grad[i]->SetDefaultUnit("tesla");
      
      
    //random limit	
    //ptot
    sprintf(strCmd,"/mydet/particle%d/momentumHigh",i+1);
    sprintf(strGuid,"Upper limit of total momentum for primary particle %d, triggered by ptot<0",i+1);
    momentumHighCmd[i] = new g4uic_double_unit(strCmd,this);
    momentumHighCmd[i]->SetGuidance(strGuid);
    momentumHighCmd[i]->SetParameterName("momentumHigh",false);
    momentumHighCmd[i]->SetDefaultUnit("GeV");	
      
    sprintf(strCmd,"/mydet/particle%d/momentumLow",i+1);
    sprintf(strGuid,"Lower limit of total momentum for primary particle %d, triggered by ptot<0",i+1);
    momentumLowCmd[i] = new g4uic_double_unit(strCmd,this);
    momentumLowCmd[i]->SetGuidance(strGuid);
    momentumLowCmd[i]->SetParameterName("momentumLow",false);
    momentumLowCmd[i]->SetDefaultUnit("GeV");
	    
    //theta
    sprintf(strCmd,"/mydet/particle%d/thetaHigh",i+1);
    sprintf(strGuid,"Upper limit of theta for primary particle %d, triggered by theta<0",i+1);
    thetaHighCmd[i] = new g4uic_double_unit(strCmd,this);
    thetaHighCmd[i]->SetGuidance(strGuid);
    thetaHighCmd[i]->SetParameterName("thetaHigh",false);
    thetaHighCmd[i]->SetRange("thetaHigh>=-10.0 && thetaHigh<=180.0");
    thetaHighCmd[i]->SetDefaultUnit("deg");	
	    
    sprintf(strCmd,"/mydet/particle%d/thetaLow",i+1);
    sprintf(strGuid,"Lower limit of theta for primary particle %d, triggered by theta<0",i+1);
    thetaLowCmd[i] = new g4uic_double_unit(strCmd,this);
    thetaLowCmd[i]->SetGuidance(strGuid);
    thetaLowCmd[i]->SetParameterName("thetaLow",false);
    thetaLowCmd[i]->SetRange("thetaLow>=-10.0 && thetaLow<=180.0");
    thetaLowCmd[i]->SetDefaultUnit("deg");
	    
    //phi
    sprintf(strCmd,"/mydet/particle%d/phiHigh",i+1);
    sprintf(strGuid,"Upper limit of phi for primary particle %d, triggered by phi<-360",i+1);
    phiHighCmd[i] = new g4uic_double_unit(strCmd,this);
    phiHighCmd[i]->SetGuidance(strGuid);
    phiHighCmd[i]->SetParameterName("phiHigh",false);
    phiHighCmd[i]->SetRange("phiHigh>=-360.0 && phiHigh<=360.0");
    phiHighCmd[i]->SetDefaultUnit("deg");
	    
    sprintf(strCmd,"/mydet/particle%d/phiLow",i+1);
    sprintf(strGuid,"Lower limit of phi for primary particle %d, triggered by phi<-360",i+1);
    phiLowCmd[i] = new g4uic_double_unit(strCmd,this);
    phiLowCmd[i]->SetGuidance(strGuid);	
    phiLowCmd[i]->SetParameterName("phiLow",false);
    phiLowCmd[i]->SetRange("phiLow>=-360.0 && phiLow<=360.0");
    phiLowCmd[i]->SetDefaultUnit("deg");
	    
    //sigma ptot, theta, phi	
    //ptot
    sprintf(strCmd,"/mydet/particle%d/sigmaMomentum",i+1);
    sprintf(strGuid,"Sigma total momentum for primary particle %d, will be used only if user specify a valid ptot(ptot>0)",i+1);
    sigmaMomCmd[i] = new g4uic_double_unit(strCmd,this);
    sigmaMomCmd[i]->SetGuidance(strGuid);
    sigmaMomCmd[i]->SetParameterName("sigmaPtot",false);
    sigmaMomCmd[i]->SetRange("sigmaPtot>=0.0");
    sigmaMomCmd[i]->SetDefaultUnit("GeV");	
    //theta
    sprintf(strCmd,"/mydet/particle%d/sigmaTheta",i+1);
    sprintf(strGuid,"Sigma theta for primary particle %d, will be used only if user specify a valid theta(0<=theta<=180)",i+1);
    sigmaThetaCmd[i] = new g4uic_double_unit(strCmd,this);
    sigmaThetaCmd[i]->SetGuidance(strGuid);
    sigmaThetaCmd[i]->SetParameterName("sigmaTheta",false);
    sigmaThetaCmd[i]->SetRange("sigmaTheta>=0.0");
    sigmaThetaCmd[i]->SetDefaultUnit("deg");	
    //phi
    sprintf(strCmd,"/mydet/particle%d/sigmaPhi",i+1);
    sprintf(strGuid,"Sigma phi for primary particle %d, will be used only if user specify a valid phi(-360<=phi<=360)",i+1);
    sigmaPhiCmd[i] = new g4uic_double_unit(strCmd,this);
    sigmaPhiCmd[i]->SetGuidance(strGuid);
    sigmaPhiCmd[i]->SetParameterName("sigmaPhi",false);
    sigmaPhiCmd[i]->SetRange("sigmaPhi>=0.0");
    sigmaPhiCmd[i]->SetDefaultUnit("deg");
	    
    // particle momentun ptot,theta,phi or p3V	
    //ptot
    sprintf(strCmd,"/mydet/particle%d/momentum",i+1);
    sprintf(strGuid,"Center total momentum for primary particle %d.",i+1);
    momentumCmd[i] = new g4uic_double_unit(strCmd,this);
    momentumCmd[i]->SetGuidance(strGuid);
    momentumCmd[i]->SetGuidance("if ptot<0 then generate ptot in [ptotlow,ptothigh],otherwise in range [ptot-sigmaptot,ptot+sigmaptot]");
    momentumCmd[i]->SetParameterName("ptot",false);
    momentumCmd[i]->SetDefaultUnit("GeV");
    //theta
    sprintf(strCmd,"/mydet/particle%d/theta",i+1);
    sprintf(strGuid,"Center theta angle for primary particle %d",i+1);
    thetaCmd[i] = new g4uic_double_unit(strCmd,this);
    thetaCmd[i]->SetGuidance(strGuid);
    thetaCmd[i]->SetGuidance("if theta<0 then generate theta in range [thetalow,thetahigh],otherwise in range [theta-sigmatheta,ptot+sigmatheta]");
    thetaCmd[i]->SetParameterName("theta",false);
    thetaCmd[i]->SetDefaultUnit("deg");
    //phi
    sprintf(strCmd,"/mydet/particle%d/phi",i+1);
    sprintf(strGuid,"Center phi for primary particle %d",i+1);
    phiCmd[i] = new g4uic_double_unit(strCmd,this);
    phiCmd[i]->SetGuidance(strGuid);
    phiCmd[i]->SetGuidance("if phi<-360 then generate phi in [philow,phihigh],otherwise in [phi-sigmaphi,phi+sigmaphi]");
    phiCmd[i]->SetParameterName("phi",false);
    phiCmd[i]->SetDefaultUnit("deg");
    //theta_tr
    sprintf(strCmd,"/mydet/particle%d/theta_tr",i+1);
    sprintf(strGuid,"Center theta_tr angle for primary particle %d",i+1);
    theta_trCmd[i] = new g4uic_double_unit(strCmd,this);
    theta_trCmd[i]->SetGuidance(strGuid);
    theta_trCmd[i]->SetGuidance("if theta_tr<0 then generate theta in range [thetalow,thetahigh],otherwise in range [theta-sigmatheta,ptot+sigmatheta]");
    theta_trCmd[i]->SetParameterName("theta_tr",false);
    theta_trCmd[i]->SetDefaultUnit("deg");
    //phi_tr
    sprintf(strCmd,"/mydet/particle%d/phi_tr",i+1);
    sprintf(strGuid,"Center phi_tr for primary particle %d",i+1);
    phi_trCmd[i] = new g4uic_double_unit(strCmd,this);
    phi_trCmd[i]->SetGuidance(strGuid);
    phi_trCmd[i]->SetGuidance("if phi_tr<-360 then generate phi in [philow,phihigh],otherwise in [phi-sigmaphi,phi+sigmaphi]");
    phi_trCmd[i]->SetParameterName("phi_tr",false);
    phi_trCmd[i]->SetDefaultUnit("deg");

    //theta_ctr
    sprintf(strCmd,"/mydet/particle%d/theta_ctr",i+1);
    sprintf(strGuid,"Center theta_ctr angle for primary particle %d",i+1);
    theta_ctrCmd[i] = new g4uic_double_unit(strCmd,this);
    theta_ctrCmd[i]->SetGuidance(strGuid);
    theta_ctrCmd[i]->SetGuidance("if theta_ctr<0 then generate theta in range [thetalow,thetahigh],otherwise in range [theta-sigmatheta,ptot+sigmatheta]");
    theta_ctrCmd[i]->SetParameterName("theta_ctr",false);
    theta_ctrCmd[i]->SetDefaultUnit("deg");
    //phi_ctr
    sprintf(strCmd,"/mydet/particle%d/phi_ctr",i+1);
    sprintf(strGuid,"Center phi_ctr for primary particle %d",i+1);
    phi_ctrCmd[i] = new g4uic_double_unit(strCmd,this);
    phi_ctrCmd[i]->SetGuidance(strGuid);
    phi_ctrCmd[i]->SetGuidance("if phi_ctr<-360 then generate phi in [philow,phihigh],otherwise in [phi-sigmaphi,phi+sigmaphi]");
    phi_ctrCmd[i]->SetParameterName("phi_ctr",false);
    phi_ctrCmd[i]->SetDefaultUnit("deg");
    //new momentun3V cmd
    sprintf(strCmd,"/mydet/particle%d/momentum3V",i+1);
    sprintf(strGuid,"Set 3-vector momentum for primary particle %d",i+1);
    momentum3VCmd[i] = new G4UIcmdWith3VectorAndUnit(strCmd,this);
    momentum3VCmd[i]->SetGuidance(strGuid);   
    momentum3VCmd[i]->SetParameterName("px","py","pz",false);
    momentum3VCmd[i]->SetDefaultUnit("GeV");

    //particle type cmds
    sprintf(strCmd,"/mydet/particle%d/particleName",i+1);
    sprintf(strGuid,"PGD particle name for primary particle %d",i+1);
    particleNameCmd[i] = new G4UIcmdWithAString(strCmd,this);
    particleNameCmd[i]->SetGuidance(strGuid); 
    particleNameCmd[i]->SetParameterName("particleName",false);

    //add PDGCode command	  
    sprintf(strCmd,"/mydet/particle%d/particlePDGCode",i+1);
    sprintf(strGuid,"Set the PDG code for primary particle %d",i+1);
    particlePDGCodeCmd[i] = new g4uic_int(strCmd,this);
    particlePDGCodeCmd[i]->SetGuidance(strGuid);  	 
    particlePDGCodeCmd[i]->SetParameterName("particlePDGcode",false);

    //event generator engine type cmds
    sprintf(strCmd,"/mydet/particle%d/engine",i+1);
    sprintf(strGuid,"Specify the event generator engine for primary particle %d",i+1);
    engineTypeCmd[i] = new G4UIcmdWithAString(strCmd,this);
    engineTypeCmd[i]->SetGuidance(strGuid); 		
    engineTypeCmd[i]->SetGuidance("Candidates includes the following: Uniform, HRSElasEl, HRSElasNucleus, HRSQuasiElasEl,");
    engineTypeCmd[i]->SetGuidance("HRSQusiElasNucleon, BoNuSProton, FastProton, RootHisto, RootNtuple, H90UserFit, Compton, TowBody");
    engineTypeCmd[i]->SetGuidance("For elastic, quasi-elastic and compon engines, please set electron or photon as particle 1,");
    engineTypeCmd[i]->SetGuidance("and set the outgoing target (nucleon or nucleus) as partice 2-8.");
    engineTypeCmd[i]->SetParameterName("engineType",false);

    //to randomize in TCS or HCS		
    sprintf(strCmd,"/mydet/particle%d/detectorAngle",i+1);
    sprintf(strGuid,"Set the detector angle primary particle %d",i+1);
    detectorAngleCmd[i] = new g4uic_double_unit(strCmd,this);
    detectorAngleCmd[i]->SetGuidance(strGuid);  	 
    detectorAngleCmd[i]->SetGuidance("This Angle will be used if randomizeInTCS is true");  	 
    detectorAngleCmd[i]->SetParameterName("detectorAngle",false);
    detectorAngleCmd[i]->SetDefaultUnit("deg");

    sprintf(strCmd,"/mydet/particle%d/randomizeInTCS",i+1);
    sprintf(strGuid,"Randomize angles in Transportation system or Hall system for primary particle %d",i+1);
    randmizeInTCSCmd[i] = new g4uic_int(strCmd,this);
    randmizeInTCSCmd[i]->SetGuidance(strGuid); 
    randmizeInTCSCmd[i]->SetGuidance("if randomizeInTCS<=0, then randomize in HCS,otherwise randomize in TCS."); 
    randmizeInTCSCmd[i]->SetGuidance("One can choose the radomize shape by giving a value: 1 is rectangle 2 is eplise.");
    randmizeInTCSCmd[i]->SetGuidance("By default use HCS. This comand will be ignored if randomize is not used.");
    randmizeInTCSCmd[i]->SetParameterName("randomizeInTCS",false);

    //theta_tr
    sprintf(strCmd,"/mydet/particle%d/outPlaneAngleHigh",i+1);
    sprintf(strGuid,"Upper limit of theta_tr (out-of-plane angle) in TCS for primary particle %d, triggered by theta<0",i+1);
    outPlaneAngleHighCmd[i] = new g4uic_double_unit(strCmd,this);
    outPlaneAngleHighCmd[i]->SetGuidance(strGuid);
    outPlaneAngleHighCmd[i]->SetParameterName("outPlaneAngleHigh",false);
    outPlaneAngleHighCmd[i]->SetRange("outPlaneAngleHigh>=-180.0 && outPlaneAngleHigh<=180.0");
    outPlaneAngleHighCmd[i]->SetDefaultUnit("deg");	

    sprintf(strCmd,"/mydet/particle%d/outPlaneAngleLow",i+1);
    sprintf(strGuid,"Lower limit of theta_tr (out-of-plane angle) in TCS for primary particle %d, triggered by theta<0",i+1);
    outPlaneAngleLowCmd[i] = new g4uic_double_unit(strCmd,this);
    outPlaneAngleLowCmd[i]->SetGuidance(strGuid);
    outPlaneAngleLowCmd[i]->SetParameterName("outPlaneAngleLow",false);
    outPlaneAngleLowCmd[i]->SetRange("outPlaneAngleLow>=-180.0 && outPlaneAngleLow<=180.0");
    outPlaneAngleLowCmd[i]->SetDefaultUnit("deg");
		
    //phi_tr
    sprintf(strCmd,"/mydet/particle%d/inPlaneAngleHigh",i+1);
    sprintf(strGuid,"Upper limit of phi_tr (in-plane angle) in TCS for primary particle %d, triggered by phi<-360",i+1);
    inPlaneAngleHighCmd[i] = new g4uic_double_unit(strCmd,this);
    inPlaneAngleHighCmd[i]->SetGuidance(strGuid);
    inPlaneAngleHighCmd[i]->SetParameterName("inPlaneAngleHigh",false);
    inPlaneAngleHighCmd[i]->SetRange("inPlaneAngleHigh>=-180.0 && inPlaneAngleHigh<=180.0");
    inPlaneAngleHighCmd[i]->SetDefaultUnit("deg");

    sprintf(strCmd,"/mydet/particle%d/inPlaneAngleLow",i+1);
    sprintf(strGuid,"Lower limit of phi_tr (in-plane angle) in TCS for primary particle %d, triggered by phi<-360",i+1);
    inPlaneAngleLowCmd[i] = new g4uic_double_unit(strCmd,this);
    inPlaneAngleLowCmd[i]->SetGuidance(strGuid);	
    inPlaneAngleLowCmd[i]->SetParameterName("inPlaneAngleLow",false);
    inPlaneAngleLowCmd[i]->SetRange("inPlaneAngleLow>=-180.0 && inPlaneAngleLow<=180.0");
    inPlaneAngleLowCmd[i]->SetDefaultUnit("deg");
		
    if(i>0) {
      sprintf(strCmd,"/mydet/particle%d/coupleToPrimary",i+1);
      sprintf(strGuid,"For two-body reactions, couple primary %d to another primary with given ordinal number",i+1);
      coupleToPrimaryCmd[i] = new g4uic_int(strCmd,this);
      coupleToPrimaryCmd[i]->SetGuidance(strGuid);	
      coupleToPrimaryCmd[i]->SetGuidance("The fact is to use momentum conservation formula to set Px_this=-Px_coupled");
      coupleToPrimaryCmd[i]->SetGuidance("Py_this=-Py_coupled and Pz_this=beamEnergy-Pz_coupled");
      coupleToPrimaryCmd[i]->SetGuidance("Note that the provided primary ordinal number should less than current ordinal number,");	
      coupleToPrimaryCmd[i]->SetGuidance("otherwise will not be coupled to any primary. To cancel coupling, one just need to set the");	
      coupleToPrimaryCmd[i]->SetGuidance("coupled primary ordinal number to 0 or any number larger than current ordinal number.");	
      coupleToPrimaryCmd[i]->SetParameterName("coupleToPrimary",false);
      coupleToPrimaryCmd[i]->SetRange("0 <= coupleToPrimary ");
      coupleToPrimaryCmd[i]->SetDefaultValue(0);
    }
  }

  randomCmd = new G4UIcmdWithABool("/mydet/randomizePrimary",this);
  randomCmd->SetGuidance("Boolean flag for randomizing primary particle types.");
  randomCmd->SetGuidance("In case this flag is false, you can select the primary particle by the following cmds");
  randomCmd->SetGuidance("/gun/particle; /mydet/particle#/particleName; /mydet/particle#/particlePDGCode; where # can be 1-8");
  randomCmd->SetParameterName("randomizeprimary",false);
  randomCmd->SetDefaultValue(false);

  particleNumCmd = new g4uic_int("/mydet/particleNum",this);
  particleNumCmd->SetGuidance("Number of primary particles in this event");
  sprintf(strGuid,"No more than %d primary articles.",MaxPrimaryNum);
  particleNumCmd->SetGuidance(strGuid);
  particleNumCmd->SetParameterName("particleNum",false);	
  sprintf(strGuid,"particleNum>=1 && particleNum<=%d",MaxPrimaryNum); 
  particleNumCmd->SetRange(strGuid);
  particleNumCmd->SetDefaultValue(1);

  //vertex
  char strRandomTrig[512];
  sprintf(strRandomTrig,"if(z<%f mm|| z>%f mm) && vertexFlag<3), a random z in range [zlow,zhigh] will be generated. \n \
						  vertexFlag will be set to: 3 for fixed point 3V; 2 for Random BeamLine and 1 for Ramdom RZ",
	  mpga->kZLowTrig/mm,mpga->kZHighTrig/mm);
	
  //random
  gunZLowCmd = new g4uic_double_unit("/mydet/gunZLow",this);
  gunZLowCmd->SetGuidance("Lower limit of Z position in the vertex random generator");	
  gunZLowCmd->SetGuidance(strRandomTrig);
  gunZLowCmd->SetParameterName("zlow",false);
  gunZLowCmd->SetDefaultValue(0.);
  gunZLowCmd->SetDefaultUnit("mm");
  //sprintf(strGuid,"zlow>=%f && zlow<=%f",mpga->kZLow/mm,mpga->kZHigh/mm);  //zlow>-15. && zlow<15.

  // Random generator of coordinate on sieve
  //random X (mm) - low
  sieveXLowCmd = unique_ptr<g4uic_double_unit>(new g4uic_double_unit("/mydet/sieveXLow",this));
  sieveXLowCmd->SetParameterName("sieve_x_low",false);
  sieveXLowCmd->SetDefaultValue(0.);
  sieveXLowCmd->SetDefaultUnit("mm");
  //random X (mm) - high
  sieveXHighCmd = unique_ptr<g4uic_double_unit>(new g4uic_double_unit("/mydet/sieveXHigh",this));
  sieveXHighCmd->SetParameterName("sieve_x_high",false);
  sieveXHighCmd->SetDefaultValue(0.);
  sieveXHighCmd->SetDefaultUnit("mm");

  //random Y (mm) - low
  sieveYLowCmd = unique_ptr<g4uic_double_unit>(new g4uic_double_unit("/mydet/sieveYLow",this));
  sieveYLowCmd->SetParameterName("sieve_y_low",false);
  sieveYLowCmd->SetDefaultValue(0.);
  sieveYLowCmd->SetDefaultUnit("mm");
  //random Y (mm) - high
  sieveYHighCmd = unique_ptr<g4uic_double_unit>(new g4uic_double_unit("/mydet/sieveYHigh",this));
  sieveYHighCmd->SetParameterName("sieve_y_high",false);
  sieveYHighCmd->SetDefaultValue(0.);
  sieveYHighCmd->SetDefaultUnit("mm");
  
  
  
  // X coordinate of vertex
  //random X (mm) - low
  gunXLowCmd = unique_ptr<g4uic_double_unit>(new g4uic_double_unit("/mydet/gunXLow",this));
  gunXLowCmd->SetParameterName("x_low",false);
  gunXLowCmd->SetDefaultValue(0.);
  gunXLowCmd->SetDefaultUnit("mm");
  
  //random X (mm) - high
  gunXHighCmd = unique_ptr<g4uic_double_unit>(new g4uic_double_unit("/mydet/gunXHigh",this));
  gunXHighCmd->SetParameterName("x_high",false);
  gunXHighCmd->SetDefaultValue(0.);
  gunXHighCmd->SetDefaultUnit("mm");

  // Y coordinate of vertex
  //random Y (mm) - low
  gunYLowCmd = unique_ptr<g4uic_double_unit>(new g4uic_double_unit("/mydet/gunYLow",this));
  gunYLowCmd->SetParameterName("y_low",false);
  gunYLowCmd->SetDefaultValue(0.);
  gunYLowCmd->SetDefaultUnit("mm");
  
  //random Y (mm) - high
  gunYHighCmd = unique_ptr<g4uic_double_unit>(new g4uic_double_unit("/mydet/gunYHigh",this));
  gunYHighCmd->SetParameterName("y_high",false);
  gunYHighCmd->SetDefaultValue(0.);
  gunYHighCmd->SetDefaultUnit("mm");
  
 
  //gunZLowCmd->SetRange(strGuid);

  gunZHighCmd = new g4uic_double_unit("/mydet/gunZHigh",this);
  gunZHighCmd->SetGuidance("Higher limit of Z position in the vertex random generator");
  gunZHighCmd->SetGuidance(strRandomTrig);
  gunZHighCmd->SetParameterName("zhigh",false);
  gunZHighCmd->SetDefaultValue(0.);
  gunZHighCmd->SetDefaultUnit("mm");

  gunXCmd = new g4uic_double_unit("/mydet/gunX",this);
  gunXCmd->SetGuidance("X position of particle gun");
  gunXCmd->SetGuidance(strRandomTrig);
  gunXCmd->SetParameterName("x",false);
  gunXCmd->SetDefaultValue(    0.0);
  gunXCmd->SetDefaultUnit("mm");

  gunYCmd = new g4uic_double_unit("/mydet/gunY",this);
  gunYCmd->SetGuidance("Y position of particle gun");
  gunYCmd->SetGuidance(strRandomTrig);
  gunYCmd->SetParameterName("y",false);
  gunYCmd->SetDefaultValue(    0.0);
  gunYCmd->SetDefaultUnit("mm");

  gunZCmd = new g4uic_double_unit("/mydet/gunZ",this);
  gunZCmd->SetGuidance("Z position of particle gun");
  gunZCmd->SetGuidance(strRandomTrig);
  gunZCmd->SetParameterName("z",false);
  gunZCmd->SetDefaultValue(-3200.0);
  gunZCmd->SetDefaultUnit("mm");

  gunX_trCmd = new g4uic_double_unit("/mydet/gunX_tr",this);
  gunX_trCmd->SetGuidance("X position of particle gun in transport coordinates");
  gunX_trCmd->SetGuidance(strRandomTrig);
  gunX_trCmd->SetParameterName("x_tr",false);
  gunX_trCmd->SetDefaultValue(    0.0);
  gunX_trCmd->SetDefaultUnit("mm");

  gunY_trCmd = new g4uic_double_unit("/mydet/gunY_tr",this);
  gunY_trCmd->SetGuidance("Y position of particle gun in transport coordinates");
  gunY_trCmd->SetGuidance(strRandomTrig);
  gunY_trCmd->SetParameterName("y_tr",false);
  gunY_trCmd->SetDefaultValue(    0.0);
  gunY_trCmd->SetDefaultUnit("mm");

  gunZ_trCmd = new g4uic_double_unit("/mydet/gunZ_tr",this);
  gunZ_trCmd->SetGuidance("Z position of particle gun in transport coordinates");
  gunZ_trCmd->SetGuidance(strRandomTrig);
  gunZ_trCmd->SetParameterName("z_tr",false);
  gunZ_trCmd->SetDefaultValue(-3200.0);
  gunZ_trCmd->SetDefaultUnit("mm");

  gunX_ctrCmd = new g4uic_double_unit("/mydet/gunX_ctr",this);
  gunX_ctrCmd->SetGuidance("X position of particle gun in transport coordinates at col");
  gunX_ctrCmd->SetGuidance(strRandomTrig);
  gunX_ctrCmd->SetParameterName("x_ctr",false);
  gunX_ctrCmd->SetDefaultValue(    0.0);
  gunX_ctrCmd->SetDefaultUnit("mm");

  gunY_ctrCmd = new g4uic_double_unit("/mydet/gunY_ctr",this);
  gunY_ctrCmd->SetGuidance("Y position of particle gun in transport coordinates at col");
  gunY_ctrCmd->SetGuidance(strRandomTrig);
  gunY_ctrCmd->SetParameterName("y_ctr",false);
  gunY_ctrCmd->SetDefaultValue(    0.0);
  gunY_ctrCmd->SetDefaultUnit("mm");

  gunZ_ctrCmd = new g4uic_double_unit("/mydet/gunZ_ctr",this);
  gunZ_ctrCmd->SetGuidance("Z position of particle gun in transport coordinates at col");
  gunZ_ctrCmd->SetGuidance(strRandomTrig);
  gunZ_ctrCmd->SetParameterName("z_ctr",false);
  gunZ_ctrCmd->SetDefaultValue(-3200.0);
  gunZ_ctrCmd->SetDefaultUnit("mm");

  gunRLowCmd = new g4uic_double_unit("/mydet/gunRLow",this);
  gunRLowCmd->SetGuidance("The lower limit of the particle gun. If raster mode is circle (1), it is the inner radius.");
  gunRLowCmd->SetGuidance("if raster mode is rectangle (2) or elipse (3), it is the half length in X axis.");
  gunRLowCmd->SetParameterName("rlow",false);
  //sprintf(strGuid,"rlow>=0 && rlow<=%.f",mpga->kRasterR/mm);  //"rlow>=0 && rlow<=20.0"
  gunRLowCmd->SetRange("rlow>=0 && rlow<=20.0");
  gunRLowCmd->SetDefaultValue(0.);
  gunRLowCmd->SetDefaultUnit("mm");

  gunRHighCmd = new g4uic_double_unit("/mydet/gunRHigh",this);
  gunRHighCmd->SetGuidance("The higher limit of the particle gun. If raster mode is circle (1), it is the outer radius.");
  gunRHighCmd->SetGuidance("if raster mode is rectangle (2) or elipse (3), it is the half length in Y axis.");
  gunRHighCmd->SetParameterName("rhigh",false);
  //sprintf(strGuid,"rhigh>=0 && rhigh<=%f",mpga->kRasterR/mm);  //"rhigh>=0 && rhigh<=20.0"
  gunRHighCmd->SetRange("rhigh>=0 && rhigh<=20.0");
  gunRHighCmd->SetDefaultValue(0.0);
  gunRHighCmd->SetDefaultUnit("mm");

  //beamline pivot point position x0,y0,z0
  fixedPointBL3VCmd = new G4UIcmdWith3VectorAndUnit("/mydet/fixedPointBL3V",this);
  fixedPointBL3VCmd->SetGuidance("Specify a pivot point (3-vector:X0,Y0,Z0) which the beam line go throught.");
  fixedPointBL3VCmd->SetGuidance("/mydet/fixedPointBL3V  must be used together with /mydet/slopeBL3V");
  fixedPointBL3VCmd->SetGuidance("Set slope for the beamline with a 3-vector (A, B, Z in mm), where A=dXdZ,B=dYdZ ==> X=X0+A*(Z-Z0) and Y=Y0+A*(Z-Z0) ");
  fixedPointBL3VCmd->SetGuidance("if cmd /mydet/fixedPointBL3V or /mydet/slopeBL3V used,  vertexFlag will be set to 2 aotumaticly");
  fixedPointBL3VCmd->SetGuidance("To unset it, you will have to use /mydet/vertexMode 1 command.");
  fixedPointBL3VCmd->SetParameterName("X0","Y0","Z0",false,true);
  fixedPointBL3VCmd->SetDefaultUnit("mm");
  fixedPointBL3VCmd->SetUnitCategory("Length");
  fixedPointBL3VCmd->SetUnitCandidates("microm mm cm m");

  //beamline slope a=dxdz,b=dydz,z, z must in mm unit
  slopeBL3VCmd = new G4UIcmdWith3Vector("/mydet/slopeBL3V",this);
  slopeBL3VCmd->SetGuidance("Set slope for the beamline in 3-vector way (A, B, Z_mm).where A=dXdZ,B=dYdZ ==> X=X0+A*(Z-Z0) and Y=Y0+A*(Z-Z0) ");
  slopeBL3VCmd->SetGuidance("Z_mm here is the vertex Z position in mm unit that will be used");	
  slopeBL3VCmd->SetGuidance(strRandomTrig);
  slopeBL3VCmd->SetGuidance("In case this cmd is using, /mydet/gunZ will be ignored");
  slopeBL3VCmd->SetGuidance("/mydet/slopeBL3V  must be used together with /mydet/fixedPointBL3V");
  slopeBL3VCmd->SetGuidance("if command /mydet/fixedPointBL3V or /mydet/slopeBL3V are used, vertexFlag will be set to 2 aotumaticly");
  slopeBL3VCmd->SetGuidance("To unset it, you will have to use /mydet/vertexMode 1 command.");
  slopeBL3VCmd->SetParameterName("A","B","Z",false,true);

  rasterModeCmd = new g4uic_int("/mydet/rasterMode",this);
  rasterModeCmd->SetGuidance("Setup the raster mode: 1 is circle, 2 is rectangle, 3 is elipse. By default it is 1, circle.");
  rasterModeCmd->SetGuidance("In circle mode (1), you can set up gunRLow as inner radius and gunRHigh as outer radius.");
  rasterModeCmd->SetGuidance("In rectangle mode (2) and elipse(3), use gunRLow is half-X and gunRHigh as hall-Y.");
  rasterModeCmd->SetGuidance("The command to set gunRLow and gunRHigh is /mydet/gunRLow and /mydet/gunRHigh.");
  rasterModeCmd->SetGuidance("If you want to turn off the raster, you just need to set both gunRLow and gunRHigh ZERO");
  rasterModeCmd->SetParameterName("rasterMode",true);
  rasterModeCmd->SetDefaultValue(1);

  vertexModeCmd = new g4uic_int("/mydet/vertexMode",this);
  vertexModeCmd->SetGuidance("Setup the vertex mode: 1 is ramdomize in RZ along Z axis, 2 is romdaomize RZ along beam line");
  vertexModeCmd->SetGuidance("3 is to use given 3-vector as vertex without raster. In mode 1 and 2, Raster will be always on.");
  vertexModeCmd->SetGuidance("Note that vertex mode will be automaticly set to 2 if if command /mydet/fixedPointBL3V or /mydet/slopeBL3V are used");
  vertexModeCmd->SetGuidance("Vertex mode will be automaticly set to 3 if command /mydet/position3V is used. ");
  vertexModeCmd->SetGuidance("If you want to turn off the raster, you just need to set both gunRLow and gunRHigh ZERO");
  vertexModeCmd->SetParameterName("vertexMode",true);
  vertexModeCmd->SetRange("vertexMode >= 1 && vertexMode<=3");
  vertexModeCmd->SetDefaultValue(1);

  //set fixed vertex by 3-vector (X,Y,Z)
  position3VCmd = new G4UIcmdWith3VectorAndUnit("/mydet/position3V",this);
  position3VCmd->SetGuidance("Specify the vertex by 3-vector(X,Y,Z). If this command is used, ");
  position3VCmd->SetGuidance("vertexFlag will be set to 3 aotumaticly and all other vertex command will be ignored.");
  position3VCmd->SetGuidance("To unset it, you will have to use /mydet/vertexMode 1 command.");
  position3VCmd->SetParameterName("X","Y","Z",false,true);
  position3VCmd->SetDefaultUnit("mm");
  position3VCmd->SetUnitCategory("Length");
  position3VCmd->SetUnitCandidates("microm mm cm m");

  //use fast spectator proton generator
  fastPsGenCmd = new G4UIcmdWithABool("/mydet/fastPsGen",this);
  fastPsGenCmd->SetGuidance("Boolean flag to use fast spectator proton generator: will generate particles(proton) inside the acceptance region");
  fastPsGenCmd->SetGuidance("In case this flag is true, all other momentum cmd for particle 1 will not be invoked.");
  fastPsGenCmd->SetGuidance("In case this flag is false, use /mydet/particle1/#### commands to set the momentum.");
  fastPsGenCmd->SetParameterName("fastPsGenFlag",true);
  fastPsGenCmd->SetDefaultValue(true);

  //to specify the beam energy and target for some engines
  beamEnergyCmd = new g4uic_double_unit("/mydet/beamEnergy",this);
  beamEnergyCmd->SetGuidance("Specify the beam energy. Some event generator engines will use it.");
  beamEnergyCmd->SetGuidance("For example, HRSElasEl, HRSElasNucleus, HRSQuasiElasEl, HRSQusiElasNucleon");
  beamEnergyCmd->SetGuidance("At the end of event, it is also used to calculate the cross section.");
  beamEnergyCmd->SetParameterName("beamEnergy",false);
  beamEnergyCmd->SetRange("beamEnergy >= 1.0e-9");
  beamEnergyCmd->SetDefaultUnit("GeV");

  tgMassCmd = new g4uic_double_unit("/mydet/targetMass",this);
  tgMassCmd->SetGuidance("Specify the target mass. Some event generator engines will use it.");
  tgMassCmd->SetGuidance("For example, HRSElasEl, HRSElasNucleus, HRSQuasiElasEl, HRSQusiElasNucleon");
  tgMassCmd->SetGuidance("At the end of event, it is also used to calculate the cross section.");
  tgMassCmd->SetParameterName("targetMass",false);
  tgMassCmd->SetRange("targetMass > 0");
  tgMassCmd->SetDefaultUnit("GeV");

  tgAtomicNumberCmd = new G4UIcmdWithADouble("/mydet/targetAtomicNumber",this);
  tgAtomicNumberCmd->SetGuidance("Specify the target atomic number. Some event generator engines will use it.");
  tgAtomicNumberCmd->SetGuidance("For example, HRSElasEl, HRSElasNucleus");
  tgAtomicNumberCmd->SetGuidance("At the end of event, it is also used to calculate the cross section.");
  tgAtomicNumberCmd->SetParameterName("targetAtomicNumber",false);
  tgAtomicNumberCmd->SetRange("targetAtomicNumber >= 0");

  //use a 3 vector to specify the beam energy and target for some engines
  beamNtarget3VCmd = new G4UIcmdWith3Vector("/mydet/beamNtarget",this);
  beamNtarget3VCmd->SetGuidance("Specify the beam energy in GeV,target mass in GeV and target atomic number.");
  beamNtarget3VCmd->SetGuidance("The following event generator engines will use them: HRSElasEl, HRSElasNucleus, HRSQuasiElasEl, HRSQusiElasNucleon");
  beamNtarget3VCmd->SetGuidance("At the end of event, they are also used to calculate the cross section.");
  beamNtarget3VCmd->SetGuidance("Note that one must input beam energy and target mass in GeV.");	
  beamNtarget3VCmd->SetParameterName("beamEnergy","targetMass","targetAtomicNumber",false);

  //to specify the left|right HRS Momentum 
  leftHRSMomentumCmd = new g4uic_double_unit("/mydet/leftHRSMomentum",this);
  leftHRSMomentumCmd->SetGuidance("Specify the momentum for left HRS.");
  leftHRSMomentumCmd->SetParameterName("leftHRSMomentum",false);
  leftHRSMomentumCmd->SetRange("leftHRSMomentum >= 0.0");
  leftHRSMomentumCmd->SetDefaultUnit("GeV");

  rightHRSMomentumCmd = new g4uic_double_unit("/mydet/rightHRSMomentum",this);
  rightHRSMomentumCmd->SetGuidance("Specify the momentum for right HRS.");
  rightHRSMomentumCmd->SetParameterName("rightHRSMomentum",false);
  rightHRSMomentumCmd->SetRange("rightHRSMomentum >= 0.0");
  rightHRSMomentumCmd->SetDefaultUnit("GeV");

  HMSMomentumCmd = new g4uic_double_unit("/mydet/HMSMomentum",this);
  HMSMomentumCmd->SetGuidance("Specify the momentum for HMS (in the right).");
  HMSMomentumCmd->SetParameterName("HMSMomentum",false);
  HMSMomentumCmd->SetRange("HMSMomentum >= 0.0");
  HMSMomentumCmd->SetDefaultUnit("GeV");


  PhotonEFracMinCmd = new G4UIcmdWithADouble("/mydet/ComptonEngine/PhotonEFracMin",this);
  PhotonEFracMinCmd->SetGuidance("Specify the minimum energy fraction of the bremsstrahlung photon in the Compton engine.");
  PhotonEFracMinCmd->SetGuidance("Only photon with energy between PhotonEFracMin and PhotonEFracMax will be generated.");
  PhotonEFracMinCmd->SetParameterName("PhotonEFracMin",false);
  PhotonEFracMinCmd->SetRange("PhotonEFracMin >= 0.0");

  PhotonEFracMaxCmd = new G4UIcmdWithADouble("/mydet/ComptonEngine/PhotonEFracMax",this);
  PhotonEFracMaxCmd->SetGuidance("Specify the maximum energy fraction of the bremsstrahlung photon in the Compton engine.");
  PhotonEFracMaxCmd->SetGuidance("Only photon with energy between PhotonEFracMax and PhotonEFracMax will be generated.");
  PhotonEFracMaxCmd->SetParameterName("PhotonEFracMax",false);
  PhotonEFracMaxCmd->SetRange("PhotonEFracMax <= 1.0");

  //G4cout<<"HRSPrimaryGeneratorMessenger() construction done!"<<G4endl;
}

HRSPrimaryGeneratorMessenger::~HRSPrimaryGeneratorMessenger()
{
  delete mydetDir;
  for(int i=0;i<MaxPrimaryNum;i++)
    {
      delete gunDir[i];
      //		delete q1grad[i];
      delete momentumHighCmd[i];
      delete momentumCmd[i];
      delete thetaLowCmd[i];
      delete thetaHighCmd[i];
      delete phiLowCmd[i];
      delete phiHighCmd[i];

      delete sigmaMomCmd[i];
      delete sigmaThetaCmd[i];
      delete sigmaPhiCmd[i];

      delete momentumLowCmd[i];
      delete thetaCmd[i];
      delete phiCmd[i];
      delete theta_trCmd[i];
      delete phi_trCmd[i];
      delete theta_ctrCmd[i];
      delete phi_ctrCmd[i];
      delete momentum3VCmd[i];

      delete particleNameCmd[i];
      delete particlePDGCodeCmd[i];

      delete engineTypeCmd[i];
      delete detectorAngleCmd[i];
      delete randmizeInTCSCmd[i];
      delete outPlaneAngleHighCmd[i];
      delete outPlaneAngleLowCmd[i];
      delete inPlaneAngleHighCmd[i];
      delete inPlaneAngleLowCmd[i];

      if(i>0) delete coupleToPrimaryCmd[i];
    }

  delete randomCmd;
  delete particleNumCmd;
  delete gunZLowCmd;
  delete gunZHighCmd;
  delete gunXCmd;
  delete gunYCmd;
  delete gunZCmd;
  delete gunX_trCmd;
  delete gunY_trCmd;
  delete gunZ_trCmd;
  delete gunX_ctrCmd;
  delete gunY_ctrCmd;
  delete gunZ_ctrCmd;
  delete gunRLowCmd;
  delete gunRHighCmd;
  delete fixedPointBL3VCmd;
  delete slopeBL3VCmd;
  delete position3VCmd;
  delete rasterModeCmd;
  delete vertexModeCmd;

  delete fastPsGenCmd;

  delete beamEnergyCmd;
  delete tgMassCmd;
  delete tgAtomicNumberCmd;
  delete beamNtarget3VCmd;
  delete leftHRSMomentumCmd;
  delete rightHRSMomentumCmd;
  delete HMSMomentumCmd;

  delete PhotonEFracMinCmd;
  delete PhotonEFracMaxCmd;

  //G4cout<<"delete HRSPrimaryGeneratorMessenger ... done!"<<G4endl;
}

void HRSPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  
  if ( command == isRHRS_Cmd.get() ) 
    target->Set_isRHRS(isRHRS_Cmd->GetNewBoolValue(newValue)); 

  if ( command == outfileCmd.get() )
    target->Set_outfile_path(newValue); 
  
  for(int i=0;i<MaxPrimaryNum;i++) {
    //random ptot
    if( command==momentumLowCmd[i] )
      { target->SetMomentumLow(i,momentumLowCmd[i]->GetNewDoubleValue(newValue)); }
    if( command==momentumHighCmd[i] )
      { target->SetMomentumHigh(i,momentumHighCmd[i]->GetNewDoubleValue(newValue)); }
    
    //random theta
    if( command==thetaLowCmd[i] )
      { target->SetThetaLow(i,thetaLowCmd[i]->GetNewDoubleValue(newValue)); }
    if( command==thetaHighCmd[i] )
      { target->SetThetaHigh(i,thetaHighCmd[i]->GetNewDoubleValue(newValue)); }
    
    //		if( command==q1grad[i] )
    //		{ cout<<"generator q1="<<q1grad[i]->GetNewDoubleValue(newValue)/tesla<<endl; }
    //		if( command==thetaLowCmd[i] )
    //		{ cout<<"generator q1="<<thetaLowCmd[i]<<" or "<<thetaHighCmd[i]->GetNewDoubleValue(newValue)<<endl; }
    
    
    
    //random phi
    if( command==phiLowCmd[i] )
      { target->SetPhiLow(i,phiLowCmd[i]->GetNewDoubleValue(newValue)); }
    if( command==phiHighCmd[i] )
      { target->SetPhiHigh(i,phiHighCmd[i]->GetNewDoubleValue(newValue)); }
    
    //sigma
    if( command==sigmaMomCmd[i] )
      { target->SetSigmaMomentum(i,sigmaMomCmd[i]->GetNewDoubleValue(newValue)); }
    if( command==sigmaThetaCmd[i] )
      { target->SetSigmaTheta(i,sigmaThetaCmd[i]->GetNewDoubleValue(newValue)); }
    if( command==sigmaPhiCmd[i] )
      { target->SetSigmaPhi(i,sigmaPhiCmd[i]->GetNewDoubleValue(newValue)); }

    //fix p, theta, phi
    if( command==momentumCmd[i] )
      { target->SetTotalMomentum(i,momentumCmd[i]->GetNewDoubleValue(newValue)); }
    if( command==thetaCmd[i] )
      {
	target->SetTheta(i,thetaCmd[i]->GetNewDoubleValue(newValue));
	//cout << "Generating: " << scientific << setprecision(15) << thetaCmd[i]->GetNewDoubleValue(newValue) << endl;
      }
    if( command==phiCmd[i] )
      {
	target->SetPhi(i,phiCmd[i]->GetNewDoubleValue(newValue));
	//cout << "Generating: " << scientific << setprecision(15) << phiCmd[i]->GetNewDoubleValue(newValue) << endl;
      }
    
    if( command==theta_trCmd[i] ){
      G4double psi   = 5.0 * deg;
      G4double chi   =-theta_trCmd[i]->GetNewDoubleValue(newValue);
      G4double theta = acos(   cos( chi ) * cos( psi ) );
      G4double phi   = atan( - tan( chi ) / sin( psi ) );
      target->SetTheta(i,theta);
      target->SetPhi  (i,phi + pi);
      //cout << "Generating: " << scientific << setprecision(15) << theta << " " << scientific << setprecision(15) << phi + pi << endl;
    }
    if( command==phi_trCmd[i] ){
      
      target->SetTheta(i,-phi_trCmd[i]->GetNewDoubleValue(newValue) + 5.0 * pi / 180.);
      target->SetPhi  (i,pi);
      //cout << "Generating: " << scientific << setprecision(15) << phi_trCmd[i]->GetNewDoubleValue(newValue) + 5.0 * pi / 180. << " " << scientific << setprecision(15) << pi << endl;
    }
    
    if( command==theta_ctrCmd[i] ){
      G4double PivotToColFace = 1.38 * m;
      G4double psi = 12.5 * pi / 180.;
      target->SetGunX( - PivotToColFace * sin( psi ) );
      target->SetGunY( 0. );
      target->SetGunZ( + PivotToColFace * cos( psi ) );
      
      //G4double psi   = 12.5 * deg;
      G4double chi   = -theta_ctrCmd[i]->GetNewDoubleValue(newValue);
      G4double theta = acos(   cos( chi ) * cos( psi ) );
      G4double phi   = atan( - tan( chi ) / sin( psi ) );
      target->SetTheta(i,theta);
      target->SetPhi  (i,phi + pi);
      //cout << "Generating: " << scientific << setprecision(15) << theta << " " << scientific << setprecision(15) << phi + pi << endl;
    }
    if( command==phi_ctrCmd[i] ){
      G4double PivotToColFace = 1.38 * m;
      G4double psi = 12.5 * pi / 180.;
      target->SetGunX( - PivotToColFace * sin( psi ) );
      target->SetGunY( 0. );
      target->SetGunZ( + PivotToColFace * cos( psi ) );

      target->SetTheta(i,-phi_ctrCmd[i]->GetNewDoubleValue(newValue) + 12.5 * pi / 180.);
      target->SetPhi  (i,pi);
      //cout << "Generating: " << scientific << setprecision(15) << phi_ctrCmd[i]->GetNewDoubleValue(newValue) + 12.5 * pi / 180. << " " << scientific << setprecision(15) << pi << endl;
    }

    if( command==momentum3VCmd[i] )
      { target->SetMomentum3V(i,momentum3VCmd[i]->GetNew3VectorValue(newValue)); }

    //particle type
    if( command==particleNameCmd[i] )
      { target->SetParticleDefinition(i,newValue); } 
    if( command==particlePDGCodeCmd[i] )
      { target->SetParticleDefinition(i,particlePDGCodeCmd[i]->GetNewIntValue(newValue)); }

    //engine type
    if( command==engineTypeCmd[i] )
      { target->SetPrimaryEngine(i,newValue); } 

    //detector angle and randomize type
    if( command==detectorAngleCmd[i] )
      { target->SetDetectorAngle(i,detectorAngleCmd[i]->GetNewDoubleValue(newValue)); }
    if( command==randmizeInTCSCmd[i] )
      { target->SetRandomizeInTCS(i,randmizeInTCSCmd[i]->GetNewIntValue(newValue)); }
    //theta_tr
    if( command==outPlaneAngleHighCmd[i] )
      { target->SetOutPlaneAngleHigh(i,outPlaneAngleHighCmd[i]->GetNewDoubleValue(newValue)); }
    if( command==outPlaneAngleLowCmd[i] )
      { target->SetOutPlaneAngleLow(i,outPlaneAngleLowCmd[i]->GetNewDoubleValue(newValue)); }
    //phi_tr
    if( command==inPlaneAngleHighCmd[i] )
      { target->SetInPlaneAngleHigh(i,inPlaneAngleHighCmd[i]->GetNewDoubleValue(newValue)); }
    if( command==inPlaneAngleLowCmd[i] )
      { target->SetInPlaneAngleLow(i,inPlaneAngleLowCmd[i]->GetNewDoubleValue(newValue)); }
    //coupleToPrimary
    if(i>0)
      {
	if( command==coupleToPrimaryCmd[i] )
	  { target->SetCoupleToPrimary(i,coupleToPrimaryCmd[i]->GetNewIntValue(newValue)); }
      }
  }
  if( command==randomCmd )
    { target->SetRandomize(randomCmd->GetNewBoolValue(newValue)); }

  if( command==particleNumCmd )
    { target->SetParticleNum(particleNumCmd->GetNewIntValue(newValue)); }

  //track intercept with sieve
  //X-range
  if( command==sieveXLowCmd.get() )
    { target->SetGunSieveXLow(sieveXLowCmd->GetNewDoubleValue(newValue)); }
  if( command==sieveXHighCmd.get() )
    { target->SetGunSieveXHigh(sieveXHighCmd->GetNewDoubleValue(newValue)); }
  //Y-range
  if( command==sieveYLowCmd.get() )
    { target->SetGunSieveYLow(sieveYLowCmd->GetNewDoubleValue(newValue)); }
  if( command==sieveYHighCmd.get() )
    { target->SetGunSieveYHigh(sieveYHighCmd->GetNewDoubleValue(newValue)); }
  
  
  //vertex 
  //X-range
  if( command==gunXLowCmd.get() )
    { target->SetGunXLow(gunXLowCmd->GetNewDoubleValue(newValue)); }
  if( command==gunXHighCmd.get() )
    { target->SetGunXHigh(gunXHighCmd->GetNewDoubleValue(newValue)); }
  //Y-range
  if( command==gunYLowCmd.get() )
    { target->SetGunYLow(gunYLowCmd->GetNewDoubleValue(newValue)); }
  if( command==gunYHighCmd.get() )
    { target->SetGunYHigh(gunYHighCmd->GetNewDoubleValue(newValue)); }
  //Z-range
  if( command==gunZLowCmd )
    { target->SetGunZLow(gunZLowCmd->GetNewDoubleValue(newValue));
    }
  if( command==gunZHighCmd ) {
    target->SetGunZHigh(gunZHighCmd->GetNewDoubleValue(newValue));
  }

	
	
	
  if( command==gunXCmd )
    { target->SetGunX(gunXCmd->GetNewDoubleValue(newValue)); }
  if( command==gunYCmd )
    { target->SetGunY(gunYCmd->GetNewDoubleValue(newValue)); }
  if( command==gunZCmd )
    { target->SetGunZ(gunZCmd->GetNewDoubleValue(newValue)); }

  if( command==gunX_trCmd ){
    G4double chi = 90.0 * pi / 180.;
    G4double psi =  5.0 * pi / 180.;
    G4double myx =-gunX_trCmd->GetNewDoubleValue(newValue);
    target->SetGunX( myx * cos( psi ) * cos( chi ) );
    target->SetGunY( myx *              sin( chi ) );
    target->SetGunZ(-myx * sin( psi ) * cos( chi ) + target->GetGunZ() );
    
  }
  if( command==gunY_trCmd ){
    G4double chi = 90.0 * pi / 180.;
    G4double psi =  5.0 * pi / 180.;
    G4double myy =-gunY_trCmd->GetNewDoubleValue(newValue);
    target->SetGunX(-myy * cos( psi ) * sin( chi ) );
    target->SetGunY( myy *              cos( chi ) );
    target->SetGunZ( myy * sin( psi ) * sin( chi ) + target->GetGunZ() );
    
  }
  if( command==gunZ_trCmd ){
    //For the moment, this command does nothing.
  }
  if( command==gunX_ctrCmd ){
    G4double PivotToColFace = 1.38 * m;
    G4double chi = 90.0 * pi / 180.;
    G4double psi = 12.5 * pi / 180.;
    G4double myx = -gunX_ctrCmd->GetNewDoubleValue(newValue);
    target->SetGunX( myx * cos( psi ) * cos( chi ) - PivotToColFace * sin( psi ) );
    target->SetGunY( myx *              sin( chi ) );
    target->SetGunZ(-myx * sin( psi ) * cos( chi ) + PivotToColFace * cos( psi ) );
    //target->SetTheta(i,psi);
    //target->SetPhi  (i,pi);
    
  }
  if( command==gunY_ctrCmd ){
    G4double PivotToColFace = 1.38 * m;
    G4double chi = 90.0 * pi / 180.;
    G4double psi = 12.5 * pi / 180.;
    G4double myy = -gunY_ctrCmd->GetNewDoubleValue(newValue);
    target->SetGunX(-myy * cos( psi ) * sin( chi ) - PivotToColFace * sin( psi ) );
    target->SetGunY( myy *              cos( chi ) );
    target->SetGunZ(-myy * sin( psi ) * sin( chi ) + PivotToColFace * cos( psi ) );
    //target->SetTheta(i,psi);
    //target->SetPhi  (i,pi);
  }
  if( command==gunZ_ctrCmd ){
    //For the moment, this command does nothing.
  }

  if( command==gunRLowCmd )
    { target->SetGunRLow(gunRLowCmd->GetNewDoubleValue(newValue)); }
  if( command==gunRHighCmd )
    { target->SetGunRHigh(gunRHighCmd->GetNewDoubleValue(newValue)); }

  if( command==rasterModeCmd )
    { target->SetRasterMode(rasterModeCmd->GetNewIntValue(newValue)); }

  if( command==vertexModeCmd )
    { target->SetVertexMode(vertexModeCmd->GetNewIntValue(newValue)); }

  if( command==fixedPointBL3VCmd )
    { target->SetFixedPointBL3V(fixedPointBL3VCmd->GetNew3VectorValue(newValue)); }
  if( command==slopeBL3VCmd )
    { target->SetSlopeBL3V(slopeBL3VCmd->GetNew3VectorValue(newValue)); }

  if( command==position3VCmd )
    { target->SetPosition3V(position3VCmd->GetNew3VectorValue(newValue)); }

  //fast ps gen cmd	
  if( command==fastPsGenCmd )
    { target->SetFastPsGenFlag(fastPsGenCmd->GetNewBoolValue(newValue)); }

  //beam and target
  if( command==beamEnergyCmd )
    { target->SetBeamEnergy(beamEnergyCmd->GetNewDoubleValue(newValue)); }
  if( command==tgMassCmd )
    { target->SetTargetMass(tgMassCmd->GetNewDoubleValue(newValue)); }
  if( command==tgAtomicNumberCmd )
    { target->SetTargetAtomicNumber(tgAtomicNumberCmd->GetNewDoubleValue(newValue)); }

  //3V cmd to set beam N target
  if( command==beamNtarget3VCmd )
    { target->SetBeamNTarget3V(beamNtarget3VCmd->GetNew3VectorValue(newValue)); }

  //HRS momentum
  //	if( command==leftHRSMomentumCmd )
  //	{ target->SetLeftHRSMomentum(leftHRSMomentumCmd->GetNewDoubleValue(newValue)); }
  //	if( command==rightHRSMomentumCmd )
  //	{ target->SetRightHRSMomentum(rightHRSMomentumCmd->GetNewDoubleValue(newValue)); }
  if( command==HMSMomentumCmd )
    { target->SetHMSMomentum(HMSMomentumCmd->GetNewDoubleValue(newValue)); }

  //compton energy fraction	
  if( command==PhotonEFracMinCmd )
    { target->SetPhotonEFracMin(PhotonEFracMinCmd->GetNewDoubleValue(newValue)); }
  if( command==PhotonEFracMaxCmd )
    { target->SetPhotonEFracMax(PhotonEFracMaxCmd->GetNewDoubleValue(newValue)); }
}

G4String HRSPrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand * command)
{
	G4String cv;

	for(int i=0;i<MaxPrimaryNum;i++)
	{
		//random ptot theta phi
		if( command==momentumLowCmd[i] )
		{ cv = momentumLowCmd[i]->ConvertToString(target->GetMomentumLow(i),"GeV"); }
		if( command==momentumHighCmd[i] )
		{ cv = momentumHighCmd[i]->ConvertToString(target->GetMomentumHigh(i),"GeV"); }
		if( command==momentumCmd[i] )
		{ cv = momentumCmd[i]->ConvertToString(target->GetTotalMomentum(i),"GeV"); }

		if( command==thetaLowCmd[i] )
		{ cv = thetaLowCmd[i]->ConvertToString(target->GetThetaLow(i),"deg"); }
		if( command==thetaHighCmd[i] )
		{ cv = thetaHighCmd[i]->ConvertToString(target->GetThetaHigh(i),"deg"); }
		if( command==thetaCmd[i] )
		{ cv = thetaCmd[i]->ConvertToString(target->GetTheta(i),"deg"); }

		if( command==phiLowCmd[i] )
		{ cv = phiLowCmd[i]->ConvertToString(target->GetPhiLow(i),"deg"); }
		if( command==phiHighCmd[i] )
		{ cv = phiHighCmd[i]->ConvertToString(target->GetPhiHigh(i),"deg"); }
		if( command==phiCmd[i] )
		{ cv = phiCmd[i]->ConvertToString(target->GetPhi(i),"deg"); }

		//sigma ptot theta phi
		if( command==sigmaMomCmd[i] )
		{ cv = sigmaMomCmd[i]->ConvertToString(target->GetSigmaMomentum(i),"GeV"); }
		if( command==sigmaThetaCmd[i] )
		{ cv = sigmaThetaCmd[i]->ConvertToString(target->GetSigmaTheta(i),"deg"); }
		if( command==sigmaPhiCmd[i] )
		{ cv = sigmaPhiCmd[i]->ConvertToString(target->GetSigmaPhi(i),"deg"); }

		//fixed momentum
		if( command==momentum3VCmd[i] )
		{ cv = momentum3VCmd[i]->ConvertToString(target->GetMomentum3V(i),"GeV"); }

		//particle type
		if( command==particleNameCmd[i] )
		{ cv = target->GetParticleName(i); }	  
		if( command==particlePDGCodeCmd[i] )
		{ cv = particlePDGCodeCmd[i]->ConvertToString(target->GetParticlePDGCode(i)); }

		//engine type, 
		if( command==engineTypeCmd[i] )
		{ cv = target->GetPrimaryEngine(i); }

		//detector angle and randomize in TCS or not
		if( command==detectorAngleCmd[i] )
		{ cv = command->ConvertToString(target->GetDetectorAngle(i)); }
		if( command==randmizeInTCSCmd[i] )
		{ cv = command->ConvertToString(target->GetRandomizeInTCS(i)); }
		
		//theta_tr		
		if( command==outPlaneAngleHighCmd[i] )
		{ cv = command->ConvertToString(target->GetOutPlaneAngleHigh(i)); }
		if( command==outPlaneAngleLowCmd[i] )
		{ cv = command->ConvertToString(target->GetOutPlaneAngleLow(i)); }
		//phi_tr		
		if( command==inPlaneAngleHighCmd[i] )
		{ cv = command->ConvertToString(target->GetInPlaneAngleHigh(i)); }
		if( command==inPlaneAngleLowCmd[i] )
		{ cv = command->ConvertToString(target->GetInPlaneAngleLow(i)); }

		if(i>0)
		{
			if( command==coupleToPrimaryCmd[i] )
			{ cv = command->ConvertToString(target->GetCoupleToPrimary(i)); }
		}
	}

	if( command==randomCmd )
	{ cv = randomCmd->ConvertToString(target->GetRandomize()); }

	if( command==particleNumCmd )
	{ cv = particleNumCmd->ConvertToString(target->GetParticleNum()); }

	//vertex
	if( command==gunZLowCmd )
	{ cv = gunZLowCmd->ConvertToString(target->GetGunZLow(),"mm"); }
	if( command==gunZHighCmd )
	{ cv = gunZHighCmd->ConvertToString(target->GetGunZHigh(),"mm"); }
	if( command==gunXCmd )
	{ cv = gunXCmd->ConvertToString(target->GetGunX(),"mm"); }
	if( command==gunYCmd )
	{ cv = gunYCmd->ConvertToString(target->GetGunY(),"mm"); }
	if( command==gunZCmd )
	{ cv = gunZCmd->ConvertToString(target->GetGunZ(),"mm"); }



	if( command==gunRLowCmd )
	{ cv = gunRLowCmd->ConvertToString(target->GetGunRLow(),"mm"); }
	if( command==gunRHighCmd )
	{ cv = gunRHighCmd->ConvertToString(target->GetGunRHigh(),"mm"); }

	if( command==rasterModeCmd )
	{ cv = rasterModeCmd->ConvertToString(target->GetRasterMode()); }

	if( command==vertexModeCmd )
	{ cv = vertexModeCmd->ConvertToString(target->GetVertexMode()); }

	//Beam Line 3 vectors
	if( command==fixedPointBL3VCmd )
	{ cv = fixedPointBL3VCmd->ConvertToString(target->GetFixedPointBL3V()); }
	if( command==slopeBL3VCmd )
	{ cv = slopeBL3VCmd->ConvertToString(target->GetSlopeBL3V()); }
	if( command==position3VCmd )
	{ cv = position3VCmd->ConvertToString(target->GetPosition3V()); }

	if( command==fastPsGenCmd )
	{ cv = fastPsGenCmd->ConvertToString(target->GetFastPsGenFlag());}

	//beam and target
	if( command==beamEnergyCmd )
	{ cv=command->ConvertToString(target->GetBeamEnergy()); }
	if( command==tgMassCmd )
	{ cv=command->ConvertToString(target->GetTargetMass()); }
	if( command==tgAtomicNumberCmd )
	{ cv=command->ConvertToString(target->GetTargetAtomicNumber()); }

	//3V cmd to set beam N target, No get for this cmd
	if( command==beamNtarget3VCmd )
	{ cv=command->ConvertToString( target->GetBeamNTarget3V());}

	//HRS momentum
	if( command==leftHRSMomentumCmd )
	{ cv=command->ConvertToString(target->GetLeftHRSMomentum()); }
	if( command==rightHRSMomentumCmd )
	{ cv=command->ConvertToString(target->GetRightHRSMomentum()); }

	if( command==HMSMomentumCmd )
	{ cv=command->ConvertToString(target->GetHMSMomentum()); }

	//compton energy fraction	
	if( command==PhotonEFracMinCmd )
	{ cv=command->ConvertToString(target->GetPhotonEFracMin()); }
	if( command==PhotonEFracMaxCmd )
	{ cv=command->ConvertToString(target->GetPhotonEFracMax()); }

	return cv;
}

