// ********************************************************************
//
// $Id: HRSSteppingAction.cc,v3.1 2008/3/16 HRS Exp $
//
//..............................................................................
#include <iomanip>
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4String.hh"
#include "HRSSteppingAction.hh"
#include "G4FieldManager.hh"
#include "G4SteppingManager.hh"
#include "HRSSteppingActionMessenger.hh"
#include "HRSTrackInformation.hh"
#include "UsageManager.hh"
#include "HRSTransform_TCSNHCS.hh"
#include "HRSEMField.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "HRSPrimaryGeneratorAction.hh"
#include <memory>
#include "HRSCoordinate.hh"
#include "TrackData_t.hh" 

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "TROOT.h"
#include "TVector3.h"

#include "G4Electron.hh"
#include "G4Positron.hh"

#include "ApexTargetGeometry.hh"

#include <stdio.h>

extern G4VPhysicalVolume *physVol_Q1_left;
extern G4VPhysicalVolume *physVol_Q1_right;
//extern G4VPhysicalVolume *physVol_Sieve_left;
//extern G4VPhysicalVolume *physVol_Sieve_right;

extern TFileHandler *fOutFile; 

using namespace std; 
//..............................................................................

//#define STEPPING_DEBUG 2
extern UsageManager* gConfig;
//extern bool CreateFileName(char *instr,char *outstr,int type);
HRSSteppingAction::HRSSteppingAction()
{
  
  
  PrintHeadFlag=true;
  //verboseLevel=2;
  verboseLevel=0;
  vPrintList.push_back("all"); 

  
  //open the TFile. this parameter is set by '/stepping/outfile'. see the associated messenger class. 
  //fOutFile = unique_ptr<TFileHandler>(new TFileHandler(fOutFile_path)); 
  
  
  int pBookTrees=0, pBookHistos=0, pBookTxt=0;
  gConfig->GetParameter("BookTrees",pBookTrees);
  gConfig->GetParameter("SeptumOn",mSeptumOn);
  gConfig->GetParameter("LHRSAngle",mLHRSAngle);
  mLHRSAngle*=deg;
  //    mLHRSAngle=45.*deg;
  //        ang_fp=45.032*deg;
  gConfig->GetParameter("LFocalPlaneAngle",ang_fp);
  ang_fp*=deg;
  gConfig->GetParameter("BookHistos",pBookHistos);
  gConfig->GetParameter("BookTxt",pBookTxt);
  iNoSecondary=1;
  gConfig->GetArgument("NoSecondary",iNoSecondary);

  
  //fTrackData = unique_ptr<TrackData_t>(new TrackData_t); 
  
  CreateTxt=(pBookTxt>0);
  CreateRootNt=(pBookTrees>0 || pBookHistos>0);

  messenger = new HRSSteppingActionMessenger(this);//unique_ptr<HRSSteppingActionMessenger>(new HRSSteppingActionMessenger(this));
  
  //G4cout<<"HRSSteppingAction() construction done!"<<G4endl;
  prevevt_id=-1;

  //G4cout<<"HRSSteppingAction() construction done!"<<G4endl;

  //        output_vdc.open("vdc_coord.C", ios::out | ios::app );
  //        output_vdc <<" {"<<endl;
  //        output_vdc <<" TH2D * vdc = new TH2D(\"vdc xy\", \"vdc xy\", 400, -40., 40., 400, -40., 40.);"<<endl;
  
  //output_focpl.open("FocalPlanecoords.txt", ios::out | ios::app );
  
  i_st=0;
  
  //	cout<<"stepping test 1"<<endl;
  //        G4FieldManager *theFieldManager;// = theStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetFieldManager();
  //	cout<<"stepping test 1"<<endl;
  //	if(!theFieldManager) theFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  //	if(G4TransportationManager::GetTransportationManager()->GetFieldManager()->DoesFieldExist());
  //	double Point[4]={700.*cm,0.*cm,700.*cm,0.};
  //	double fieldArr[6]={0,0,0,0,0,0};
  
  P1x=0.;
  P1y=0.;
  P1z=0.;
  P2x=0.;
  P2y=0.;
  P2z=0.;
  
  //cout << "make txt output : " << (Do_TXToutput()?"true":"false") << endl; 
  
  /*
    G4FieldManager *globalFieldManager;// = theStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetFieldManager();
    G4TransportationManager *transportMgr= G4TransportationManager::GetTransportationManager();
    globalFieldManager = transportMgr->GetFieldManager();
    G4double minEps= 1.0e-5;  //   Minimum & value for smallest steps
    G4double maxEps= 1.0e-4;  //   Maximum & value for largest steps
    globalFieldManager->SetMinimumEpsilonStep( minEps );  
    globalFieldManager->SetMaximumEpsilonStep( maxEps );
    globalFieldManager->SetDeltaOneStep( 0.5e-3 * mm );  // 0.5 micrometer
    G4cout << "EpsilonStep: set min= " << minEps << " max= " << maxEps << G4endl;
  */
}

//..............................................................................

HRSSteppingAction::~HRSSteppingAction()
{
  if (GetVerbose()>=2) cout << "Calling <HRSSteppingAction::~()>..." << endl; 

  //fOutFile->CloseFile();
  //delete fOutFile; 
  
//        output_vdc <<" vdc->Draw(\"colz\");"<<endl;
//        output_vdc <<" }"<<endl;
//        output_vdc <<endl;
//        output_vdc.close();
  
  if (CreateTxt)
    {
      OutTxt.close();
    }
  vPrintList.clear();
  //G4cout<<"delete HRSSteppingAction ... done!"<<G4endl;

  if (GetVerbose()>=2) cout << "<HRSSteppingAction>: Deleting TrackData_t..." << flush;
  //delete fTrackData;
  if (GetVerbose()>=2) cout << "done!" << endl;
  
  
  if (GetVerbose()>=2) cout << "<HRSSteppingAction>: Deleting messenger..." << flush;
  delete messenger;
  if (GetVerbose()>=2) cout << "done!" << endl;
  
  if (GetVerbose()>=2) cout << "ending call of <HRSSteppingAction::~()>" << endl; 
}
//..............................................................................

void HRSSteppingAction::InitOutTxt()
{
  if (!CreateTxt) return;

  char strFileName[200], strRawName[200];
	
  std::string pOutFileName=gConfig->GetArgument("OutFileName");
  /*
    if(gHRSTree)
    {
    pOutFileName=gHRSTree->GetRootFileName();
    }
    if(pOutFileName.length()>3) 
    {
    string tmpStr(pOutFileName);
    size_t pos=tmpStr.rfind(".root");
    if(pos!=string::npos) tmpStr.replace(pos,5,".txt"); 
    strcpy(strFileName,tmpStr.c_str());
    }
    else
    {
    //sprintf(strRawName,"G4Sim_txt_%02d.txt",gHRSTree->iRunNumber);
    //CreateFileName(strRawName,strFileName);
    //in order to make the txt output has the same name with the root output, I check 
    //the existance of the ntuple root file instead of txt file
    sprintf(strRawName,"G4Sim_nt_%02d.root",gHRSTree->iRunNumber);
    CreateFileName(strRawName,strFileName,2);
    }
    OutTxt.open (strFileName, fstream::out | fstream::app);
  */
}

void HRSSteppingAction::UserSteppingAction(const G4Step* theStep)
{
  //HRSCoordinate::Arm arm = HRSCoordinate::kLeft;

  //a simple helper function to convert a G4ThreeVector to a TVector3.
  // This is because TVector3's can automatically be included in a TTree as a branch,
  // but the same is not true for G4ThreeVector-s. There is almost certainly a
  // more elegant way to do this, which probably involves making the 'G4ThreeVector'
  // class compatable with TTree-s, but I don't really care to do all that extra effort.
  auto Get_TVector3_from_G4ThreeVector = [](TVector3 &out, G4ThreeVector in)
  {
    out.SetXYZ( in.x(), in.y(), in.z() ); 
  };
  
    
  //the track we will be measuring / operating on 
  G4Track *track = theStep->GetTrack(); 
  
    
  //this struct keeps information about the track that we will need
  TrackData_t *track_data = fOutFile->Get_TData(); 

  
  
  
  //this helper function makes sure that all values are reset
  auto kill_track = [&track, this](bool write_track=false)
  {
    track->SetTrackStatus(fStopAndKill);
    //this->fOutFile->ClearTrack(write_track); 
  };

  
  //get the event number. check if this is a new track or not. 
  int event_id = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  
  //track_data->event_id = event_id;
  //track_data
  
  int track_id = track->GetTrackID(); 

  //check for secondary tracks
  if (iNoSecondary && (track_id != 1)) {
    kill_track();
    if (GetVerbose()>=3) cout << "secondary track killed" << endl;
    return;
  }

  
  //reset the track data, if its a new event/track. 
  if (event_id != track_data->event_id ||
      track_id != track_data->track_id) {
    
    if (GetVerbose()>=3)
      cout << "new evt/track. id_step="
	   << event_id << " id_trk=" << track_data->event_id << endl;
    
    track_data->event_id = event_id;
    track_data->track_id = track_id;
  }//*/
    
  
  //different helper function, extracts data from the track that we need.
  //this takes a ptr to a HRSCoordinate object, and fills it with the current track data.
  /*auto fill_coordinate = [&track, arm, this](HRSCoordinate *coord,
					     HRSCoordinate::Coord coord_sys = HRSCoordinate::kHCS)
  {
    coord->position = track->GetPosition();
    coord->momentum = track->GetMomentum(); 
    
    if (coord_sys == HRSCoordinate::kTCS) {;} //handle TCS case
    };*/ 
  
  
  //check if this is a secondary track. If so,
  
  if (iNoSecondary && track->GetTrackID()!=1) {
    
    kill_track();
    if (GetVerbose()>=3) cout << "-(secondary track killed)" << endl;
    return;
  }
    //            cout<<"the track:"<<theStep->GetTrack()->GetTrackID()<<" is killed - secondary particle"<<endl;
    //theStep->GetTrack()->SetTrackStatus(fStopAndKill);  
  
    
    //  cout << "WE ARE IN THE USER STEPPING ACTION" << endl;
  G4VPhysicalVolume* physVol_current  = theStep->GetPreStepPoint()->GetPhysicalVolume();
  /*G4String name = physVol_current->GetName();
  name.assign(name,0,70);
  G4String name1 = physVol_current->GetName();
  name1.assign(name1,0,40);
  G4String name2 = physVol_current->GetName();
  name2.assign(name2,0,7);*/ 

  
  
  
  double edeposit   =  theStep->GetTotalEnergyDeposit();
  //e_sum+=edeposit/1000.;
  //        G4int eID = 0;
  //        const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
  //        if(evt) eID = evt->GetEventID();
  
  
  
  G4int evtNb = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  evt_id=evtNb;
  if (prevevt_id != evt_id) {
    track_id_rad_tail_foil=-1*evt_id;
    px_rad_tail_foil=-999.1;
    py_rad_tail_foil=-999.2;
    pz_rad_tail_foil=-999.3 ;
    track_id_sieve_back==-1*evt_id;
    px_sieve_back=-999.1;
    py_sieve_back=-999.2;
    pz_sieve_back=-999.3;
    x_sieve_back=-1999.1;
    y_sieve_back=-1999.2;
    z_sieve_back=-1999.3;
    
    write_trig=true;
    prevevt_id =evt_id;
  }
  
  G4double z=theStep->GetTrack()->GetPosition().z()/cm;
  G4double x=theStep->GetTrack()->GetPosition().x()/cm;
  G4double y=theStep->GetTrack()->GetPosition().y()/cm;


  double y_TCS=y;
  double thet_tcs=acos(z/sqrt(x*x+z*z))-mLHRSAngle;
  //    cout<<"thet_tcs="<<thet_tcs/deg<<endl;
  double x_TCS=sqrt(x*x+z*z)*sin(thet_tcs);
  double z_TCS=sqrt(x*x+z*z)*cos(thet_tcs);
  //        cout<<"Target CS:"<<"  angle="<<thet_tcs<<"     x="<<x_TCS<<"   y="<<y_TCS<<"   z="<<z_TCS<<endl;
  /*
  G4double x_vert=theStep->GetTrack()->GetVertexPosition().x()/cm;
  G4double y_vert=theStep->GetTrack()->GetVertexPosition().y()/cm;
  G4double z_vert=theStep->GetTrack()->GetVertexPosition().z()/cm;

  float r=sqrt(x*x + y*y);
  float ang_phi=acos(x/r) * y/fabs(y);
  float ang_theta = atan(r/z);
  float eta_tmp = -1*log(fabs(tan(ang_theta/2)))*z/fabs(z);
  double xcoord=theStep->GetPreStepPoint()->GetPosition().x();
  double ycoord=theStep->GetPreStepPoint()->GetPosition().y();
  double zcoord=theStep->GetPreStepPoint()->GetPosition().z();
  double px_step=theStep->GetPreStepPoint()->GetMomentum().x()/GeV;
  double py_step=theStep->GetPreStepPoint()->GetMomentum().y()/GeV;
  double pz_step=theStep->GetPreStepPoint()->GetMomentum().z()/GeV;
  double pt=sqrt(theStep->GetTrack()->GetMomentum().x()*theStep->GetTrack()->GetMomentum().x() + theStep->GetTrack()->GetMomentum().y()*theStep->GetTrack()->GetMomentum().y());

  double p_HCS=sqrt(px_step*px_step + py_step*py_step + pz_step*pz_step);

  double p_thet=asin(px_step/sqrt(px_step*px_step + pz_step*pz_step))-mLHRSAngle;
  double p_x_TCS=sqrt(px_step*px_step + pz_step*pz_step)*sin(p_thet);
  double p_y_TCS=py_step;
  double p_z_TCS=sqrt(px_step*px_step + pz_step*pz_step)*cos(p_thet);
  double p_phi =asin(py_step/sqrt(p_HCS*p_HCS-p_x_TCS*p_x_TCS));
        
  double px_dir_vert=theStep->GetTrack()->GetVertexMomentumDirection().x();
  double py_dir_vert=theStep->GetTrack()->GetVertexMomentumDirection().y();
  double pz_dir_vert=theStep->GetTrack()->GetVertexMomentumDirection().z();
  double thet_vert = asin(px_dir_vert/sqrt(px_dir_vert*px_dir_vert + pz_dir_vert*pz_dir_vert))-mLHRSAngle;
  double p_z_vert_tcs=sqrt(px_dir_vert*px_dir_vert + pz_dir_vert*pz_dir_vert)*cos(thet_vert);
  double phi_vert = asin(py_dir_vert/p_z_vert_tcs);



  double t_targ=sqrt(x*x+z*z)*sin(acos(z/sqrt(x*x+z*z))-mLHRSAngle);
  //	output1 <<"x["<<i_st<<"]="<<sqrt(x*x+z*z)<<";  y["<<i_st<<"]="<<t_targ<<";"<<endl;
  //	output <<"x["<<i_st<<"]="<<sqrt(x*x+z*z)<<";  y["<<i_st<<"]="<<y<<";"<<endl;
  i_st++;
  
  double z00=9.961*m;
  double y00=8.4*m;
  double R_new= sqrt( (y00-yyy)*(y00-yyy) + (zzz - z00) * (zzz - z00) );
  */ 
  
  if ( (z_TCS>620.) && (z_TCS<630.) && (sqrt(x_TCS*x_TCS + y_TCS*y_TCS)>31.) ) {
    
    kill_track(); 
    if (GetVerbose()>=2) cout<<"-(the track is killed Q2 ex)"<<flush;
    return;
    //theStep->GetTrack()->SetTrackStatus(fStopAndKill);
  }
  
  double pQ3Cen_y_tmp =  358.637  + 182./sqrt(2.)/2.;
  double pQ3Cen_z_tmp = 1702.67042  + 182./sqrt(2.)/2.;

  if ((z_TCS>-10.) &&
      (z_TCS<160.) &&
      ((fabs(x_TCS)>20.) || (fabs(y_TCS)>20.))) {
    
    if (GetVerbose()>=2) cout<<"-(the track is killed before Q1)"<<flush;
    kill_track();//theStep->GetTrack()->SetTrackStatus(fStopAndKill);
    return;
  }

  if ((z_TCS>pQ3Cen_z_tmp) &&
      (z_TCS<pQ3Cen_z_tmp+10.) && 
      (sqrt( x_TCS*x_TCS + (y_TCS-pQ3Cen_y_tmp)*(y_TCS-pQ3Cen_y_tmp)) > 45.)) {
    
    if (GetVerbose()>=2) cout<<"-(the track is killed Q3 cen)"<<flush;
    kill_track();//theStep->GetTrack()->SetTrackStatus(fStopAndKill);
    return;
  }
  
  
  if (fabs(y)>820.) {
    
    if (GetVerbose()>=2) cout<<"-(the track:"<<track->GetTrackID()<<" is killed y > 820)"<<flush;
    kill_track();//theStep->GetTrack()->SetTrackStatus(fStopAndKill);
    return;
  }
  

  //septum vacuuum box
  if (mSeptumOn &&
      ((z>-109.) && z<(160.))) {

    double inch           = 2.54;
    
    double zlength        = 173.939;
    double z_sep_cen      = zlength*tan(5.*deg)/tan(12.5*deg);
    double z_real_tar     = z_sep_cen-zlength;
    double z_sept_en_min1 = z_real_tar + (17.31*inch+22.50*inch);
    double z_sept_en_max1 = z_sept_en_min1 - 3.15*inch*tan(5.*deg);
    double x_cen_sep_en   = 3.966*inch;
    double x_width_sep_en = 3.15*inch;
    double xmin_sep_en1   = x_cen_sep_en -x_width_sep_en*0.5*cos(5.*deg);
    double xmax_sep_en1   = x_cen_sep_en +x_width_sep_en*0.5*cos(5.*deg);
    double ang_en_min_1   = 5.*deg;
    double ang_en_max_1   = 5.*deg;
    double ang_en_min_2   = 6.6*deg;
    double ang_en_max_2   = 10.8*deg;
    double length_max_1   = 50.19;//* cm; // 19.76*inch
    double length_min_1   = 52.1 ;//* cm; // 19.76*inch
    double length_min_2   = 59.82;//* cm; // 23.55*inch
    double length_max_2   = 60.3 ;//* cm; // 23.74*inch

    double xmin_sep_ex1   = xmin_sep_en1 + length_min_1 * sin(ang_en_min_1);
    double xmax_sep_ex1   = xmax_sep_en1 + length_max_1 * sin(ang_en_max_1);
    double z_sept_ex_min1 = z_sept_en_min1 + length_min_1 * cos(ang_en_min_1);
    double z_sept_ex_max1 = z_sept_en_max1 + length_max_1 * cos(ang_en_max_1);

    double xmin_sep_en2   = xmin_sep_ex1;
    double xmax_sep_en2   = xmax_sep_ex1;
    double z_sept_en_min2 = z_sept_ex_min1;
    double z_sept_en_max2 = z_sept_ex_max1;

    double xmin_sep_ex2   = xmin_sep_en2 + length_min_2 * sin(ang_en_min_2);
    double xmax_sep_ex2   = xmax_sep_en2 + length_max_2 * sin(ang_en_max_2);
    double z_sept_ex_min2 = z_sept_en_min2 + length_min_2 * cos(ang_en_min_2);
    double z_sept_ex_max2 = z_sept_en_max2 + length_max_2 * cos(ang_en_max_2);

    //            cout<<"kuku  "<<z<<"     "<<x<<"        "<<xsep_tmp<<endl;
    

    /*
      if ((z>z_sept_en_min1) && (z<z_sept_ex_min1))
      {
      double xsep_tmp = xmin_sep_en1 + (xmin_sep_ex1-xmin_sep_en1)*(z-z_sept_en_min1)/(z_sept_ex_min1-z_sept_en_min1);
      if (xsep_tmp > x) theStep->GetTrack()->SetTrackStatus(fStopAndKill);
      }

      if ((z>z_sept_en_min2) && (z<z_sept_ex_min2))
      {
      double xsep_tmp = xmin_sep_en2 + (xmin_sep_ex2-xmin_sep_en2)*(z-z_sept_en_min2)/(z_sept_ex_min2-z_sept_en_min2);
      //              cout<<"kuku  "<<z<<"     "<<x<<"        "<<xsep_tmp<<endl;
      if (xsep_tmp > x) theStep->GetTrack()->SetTrackStatus(fStopAndKill);
      }

      if ((z>z_sept_en_max1) && (z<z_sept_ex_max1))
      {
      double xsep_tmp = xmax_sep_en1 + (xmax_sep_ex1-xmax_sep_en1)*(z-z_sept_en_max1)/(z_sept_ex_max1-z_sept_en_max1);
      if (xsep_tmp < x) theStep->GetTrack()->SetTrackStatus(fStopAndKill);
      }
      if ((z>z_sept_en_max2) && (z<z_sept_ex_max2))
      {
      double xsep_tmp = xmax_sep_en2 + (xmax_sep_ex2-xmax_sep_en2)*(z-z_sept_en_max2)/(z_sept_ex_max2-z_sept_en_max2);
      if (xsep_tmp < x) theStep->GetTrack()->SetTrackStatus(fStopAndKill);
      }
    */
    
    if ( (z > z_sept_en_max1) && (z < z_sept_ex_min2) ) {
      double y_en = 2.44*inch;
      double y_ex = 4.7 *inch;
      double ysep_tmp = y_en + (z-z_sept_en_max1)/(z_sept_ex_min2-z_sept_en_max1)*(y_ex-y_en);
      //              cout<<z<<"     "<<ysep_tmp<<endl;
      if (fabs(y) > ysep_tmp ) {
	if (GetVerbose()>=2) cout<<"-(the track is killed; y > ysep_tmp)"<<flush;
	//theStep->GetTrack()->SetTrackStatus(fStopAndKill);
	kill_track();
	return;
      }
	
    }
  }//if (mSeptumOn)
  
  
  
  /*G4ThreeVector position = track->GetPosition();
  printf("-(p[% 5.f,% 5.f,% 5.f])\n",
  position(0),position(1),position(2)); */ 
  
  //        if (theStep->GetTrack()->GetTrackID()==1)
  if ((*physVol_current == *physVol_Q1_left) || 
      (*physVol_current == *physVol_Q1_right)) {

    auto track_definition = track->GetDefinition(); 
    
    //check that this is either an electron or a positron
    if ( ! (track_definition == G4Electron::Electron() ||
	    track_definition == G4Positron::Positron()) ) {
      if (GetVerbose()>=2) cout << "-(Particle type neither electron nor positron.)";
      kill_track();
      return;
    }
    
    
    bool is_RHRS(track_definition == G4Positron::Positron()); 
    
    
    if (GetVerbose()>=2) cout << "-(at Q1_" << (is_RHRS?"R)":"L)") << flush; 
    /*if (GetVerbose()>=2) printf("-(pQ1[% 1.3f,% 1.3f,% 1.3f])",
      track->GetPosition().x()*100, 
      track->GetPosition().y()*100, 
      track->GetPosition().z()*100);*/
    
    double hrs_angle = ApexTargetGeometry::Get_HRS_angle(is_RHRS); 
    /*double hrs_angle = is_RHRS ? -12.50 : 12.50; 
      hrs_angle *= CLHEP::pi / 180.;*/

    double sieve_angle = ApexTargetGeometry::Get_sieve_angle(is_RHRS); 
    /*double sieve_angle = is_RHRS ? -5.372 : 5.366;
      sieve_angle *= CLHEP::pi / 180.;*/
    
    
    //fill out data about the track~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // vertex position / momentum
    
    //this is necessary, as (for some reason) you can only get the *normalized*
    // momentum at the vertex position. 
    double T_vtx = track->GetVertexKineticEnergy();
    
    G4ThreeVector momentum_vtx = track->GetVertexMomentumDirection();
    momentum_vtx *= T_vtx; 

    //get the position of the particle's vertex (Rel. to the Apex target center!) 
    G4ThreeVector position_vtx
      = track->GetVertexPosition() - ApexTargetGeometry::Get_APEX_Target_center(); 
    
    //These are saved in Hall coordinates (no need to rotate them). 
    Get_TVector3_from_G4ThreeVector( track_data->position_vtx, position_vtx/1000. );
    Get_TVector3_from_G4ThreeVector( track_data->momentum_vtx, momentum_vtx );
    
    //Now, we're going to compute their projections at the sieve-plane. 
    position_vtx.rotateY( -sieve_angle );
    momentum_vtx.rotateY( -sieve_angle );  
    
    //In TCS (which is what SIMC is expecting), the x-axis points to the floor of the
    // hall, and the y-axis points to the left (if you're facing into the opening of
    // the spectrometer.) 
        
    position_vtx.rotateZ( CLHEP::pi/2. ); 
    momentum_vtx.rotateZ( CLHEP::pi/2. ); 
    
    //These are hard-coded offsets which reflect measurements from the APEX survey.
    //Since we've already converted to meters, These values are also in meters. 
    G4ThreeVector sieve_D0 = ApexTargetGeometry::Get_sieve_pos(is_RHRS); 
    
    /*G4ThreeVector sieve_D0(  is_RHRS ?  -1.101e-3 :  -1.301e-3,
			     is_RHRS ?  -3.885e-3 :   6.672e-3,
			     is_RHRS ? 794.609e-3 : 795.766e-3 ); */ 

    
    //Project these onto the sieve plane. in sieve coordinates, the target-facing
    //opening of the central large hole is defined as the origin.
    // the x-axis is normal to the hall floor (pointing down).
    // the z-axis is coaxial with the central large hole, pointing into the septum (away
    // from target). the y-axis, then, is parallel to the hall floor, pointing away from the
    // beam (and tangent to the face of the sieve).
    double dz,dxdz,dydz; 
    
    dz   = sieve_D0.z() - position_vtx.z();
    dxdz = momentum_vtx.x()/momentum_vtx.z();
    dydz = momentum_vtx.y()/momentum_vtx.z();
    
    
    //project this track onto the plane of the sieve's front face
    position_vtx(0) += dxdz*dz;
    position_vtx(1) += dydz*dz;
    position_vtx(2) = sieve_D0.z();
    
    //shift the position vector relative to the sieve offset. 
    position_vtx = position_vtx - sieve_D0; 

    
    //These are saved in Hall coordinates (no need to rotate them). 
    Get_TVector3_from_G4ThreeVector( track_data->position_sieve, position_vtx/1000. );
    Get_TVector3_from_G4ThreeVector( track_data->momentum_sieve, momentum_vtx );

    
    /*if (GetVerbose()>=2) printf("-(pQ1[% 1.3f,% 1.3f,% 1.3f])",
				position_vtx.x()*100, 
				position_vtx.y()*100, 
				position_vtx.z()*100);*/ 
    

    // - Q1 position / momentum
    // get the momentum & positions at the Q1 plane, then 'rotate' them by the angle
    // of the respective spectrometer.
    // these units are in meters; the G4Track::GetPosition() method gives us in mm.

    //Just to reiterate - we're dividing by 1000 to conver from mm to m.
    G4ThreeVector position_Q1 = track->GetPosition();

    //Here, the units are MeV/c
    G4ThreeVector momentum_Q1 = track->GetMomentum();
    
    /*if (GetVerbose()>=2) printf("-(pQ1[% -5.3f,% -5.3f,% -5.3f])",
      position_Q1.x()*100, 
      position_Q1.y()*100, 
      position_Q1.z()*100);*/ 
    
    
    position_Q1.rotateY( -hrs_angle ); 
    momentum_Q1.rotateY( -hrs_angle ); 

    //In TCS (which is what SIMC is expecting), the x-axis points to the floor of the
    // hall, and the y-axis points to the left (if you're facing into the opening of
    // the spectrometer.) 
    
    position_Q1.rotateZ( CLHEP::pi/2. ); 
    momentum_Q1.rotateZ( CLHEP::pi/2. ); 
    
    if (GetVerbose()>=2) printf("-(pQ1[% 5.3f,% 5.3f,% 5.3f])",
				position_Q1.x()*100., 
				position_Q1.y()*100., 
				position_Q1.z()*100.);//*/   
      
    
    //project these tracks back onto the Q1 front plane
    const double Q1_front_z = 1710.95; //in meters, in rotated-HCS
    
    dz   = position_Q1.z() - Q1_front_z;
    dxdz = momentum_Q1.x()/momentum_Q1.z(); 
    dydz = momentum_Q1.y()/momentum_Q1.z(); 
    
    position_Q1(0) += -(dz*dxdz); //x
    position_Q1(1) += -(dz*dydz); //y
    position_Q1(3)  = Q1_front_z; //z

    
    //now we're ready to fill them into th TrackData_t struct
    Get_TVector3_from_G4ThreeVector( track_data->position_Q1, position_Q1/1000. );
    Get_TVector3_from_G4ThreeVector( track_data->momentum_Q1, momentum_Q1 );

    
    //now, we record this track (tell the TFileHelper object to write it to the output file)
    track_data->event_id = event_id;
    track_data->track_id = track_id; 

    
    //write this track to our TTree
    fOutFile->WriteTrack(); 

    
    //output stream (txt) 
    if (Do_TXToutput()) {
      
      /*double e_kin = theStep->GetTrack()->GetVertexKineticEnergy()+0.511;
      double px_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().x();
      double py_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().y();
      double pz_vtx_tmp=e_kin*theStep->GetTrack()->GetVertexMomentumDirection().z();*/ 
      
      if (GetVerbose()>=2) cout << " Making txt output.." << flush; 
      
      /*std::ofstream output_q1;
      output_q1.open( "Q1_front.dat", ios::out | ios::app );
    
      output_q1 //<< precision(0)
	<<evtNb <<setw(18)
	<<theStep->GetTrack()->GetTrackID() <<setw(18)
	<<z_TCS <<setw(18)
	<<p_HCS <<setw(18)
	<<p_x_TCS/p_z_TCS <<setw(18)
	<<p_y_TCS/p_z_TCS <<setw(18)
	<<x_TCS <<setw(18)
	<<y_TCS <<setw(18)
	<<px_vtx_tmp <<setw(18)
	<<py_vtx_tmp <<setw(18)
	<<pz_vtx_tmp <<setw(18)
	<<theStep->GetTrack()->GetDefinition()->GetParticleName() <<setw(18)
	<<track_id_rad_tail_foil <<setw(18)
	<<px_rad_tail_foil <<setw(18)
	<<py_rad_tail_foil <<setw(18)
	<<pz_rad_tail_foil <<setw(18)
	<<track_id_sieve_back <<setw(18)
	<<x_sieve_back <<setw(18)
	<<y_sieve_back <<setw(18)
	<<z_sieve_back <<setw(18)
	<<px_sieve_back <<setw(18)
	<<py_sieve_back <<setw(18)
	<<pz_sieve_back <<endl;
	output_q1.close();      */
      
    }//if (Do_TXToutput()) {
    

    //kill this track; if septum is on (in apex mode, we don't simulate the HRS). 
    if (mSeptumOn && Do_killAtQ1()) {
      if (GetVerbose() >= 2) cout << "-(killed at Q1)" << flush; 
      //fOutFile->WriteTrack();
      kill_track(); //theStep->GetTrack()->SetTrackStatus(fStopAndKill);
      return;
    }			
  }
  
}

//..............................................................................

void HRSSteppingAction::SetOutFilePath(const std::string &_path)
{
  cout << "<HRSSteppingAction::SetOutFilePath>: Creating TFile..." << endl;
  //fOutFile = new TFileHandler(_path.data());
  //fOutFile = unique_ptr<TFileHandler>(new TFileHandler(_path.data()));  
}
//..............................................................................
//void HRSSteppingAction::SetOutFilePath(const std::string &_name) {;}
//..............................................................................
//..............................................................................
//..............................................................................
//..............................................................................



void HRSSteppingAction::DoPrint(const G4Step* theStep)
{
}


//..............................................................................

void HRSSteppingAction::PrintHead(const G4Step* theStep,ostream& pOut)
{
}
//..............................................................................

void HRSSteppingAction::PrintStep(const G4Step* theStep,ostream& pOut)
{	
}
//..............................................................................

void HRSSteppingAction::FillRootArray(const G4Step* theStep)
{
}

