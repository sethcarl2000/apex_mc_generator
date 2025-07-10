// ********************************************************************
//
// $Id: HRSVisAttribute.cc,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//This class must be launched from HRSDetectorConstruction and should
//be instanced only once
//

#include "HRSVisAttribute.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

using namespace std;

using g4va = G4VisAttributes;

/////////////////////////////////////////////////////////////////////
HRSVisAttribute::HRSVisAttribute()
{
	ConstructVisAttribute();
}

HRSVisAttribute::~HRSVisAttribute()
{
  /*  delete Invisible_VisAtt; 
  
  delete HallVisAtt;
  delete MagneticVisAtt;

  delete WhiteVisAtt;
  delete DarkRedVisAtt;
  delete OrangeVisAtt;
  delete GrayVisAtt;
  delete YellowVisAtt;
  delete PcbGreenVisAtt;
  delete CuBrownVisAtt;
  delete DarkBlueVisAtt;
  delete VioletVisAtt;
  delete LightYellowVisAtt;

  delete BlackVisAtt; 
  delete CableVisAtt; 

  delete RedVisAtt;
  delete YellowGreenVisAtt; 
  delete LightPurpleVisAtt; 
  delete PurpleVisAtt;
  delete DarkOrangeVisAtt;
  delete SkyBlueVisAtt; 
  delete LightGreenVisAtt; 
  delete LightBlueVisAtt; 

  delete SteelVisAtt;
  delete IronVisAtt;
  delete SilverVisAtt;  
  delete LeadVisAtt;*/
}

/////////////////////////////////////////////////////////////////////
void HRSVisAttribute::ConstructVisAttribute()
{

  //invisible
  Invisible_VisAtt = unique_ptr<g4va>(new g4va); Invisible_VisAtt->SetVisibility(false); 
  
  HallVisAtt = unique_ptr<g4va>(new g4va(G4Colour(1.0,1.0,1.0))); //White
  HallVisAtt->SetVisibility(false);
  MagneticVisAtt = unique_ptr<g4va>(new g4va(G4Colour(0.9,0.9,0.9)));   // LightGray
  MagneticVisAtt->SetVisibility(false);
  
  WhiteVisAtt = unique_ptr<g4va>(new g4va(G4Colour(1.0,1.0,1.0)));	//white
  DarkRedVisAtt = unique_ptr<g4va>(new g4va(G4Colour(1.0,0.0,0.0)));//dark red, for target-wall
  OrangeVisAtt = unique_ptr<g4va>(new g4va(G4Colour(1.0,0.5,0.0)));//orange
  GrayVisAtt = unique_ptr<g4va>(new g4va(G4Colour(0.75,0.75,0.75)));	//grey
  YellowVisAtt = unique_ptr<g4va>(new g4va(G4Colour(1.0,1.0,0.0)));	//yellow, color for Ultem
  PcbGreenVisAtt = unique_ptr<g4va>(new g4va(G4Colour(0.0,1.0,0.0)));//green, color for pcb
  CuBrownVisAtt = unique_ptr<g4va>(new g4va(G4Colour(1.0,0.5,0.5)));	//color closed to Cu
  DarkBlueVisAtt = unique_ptr<g4va>(new g4va(G4Colour(0.0,0.5,1.0))); //dark blue
  VioletVisAtt = unique_ptr<g4va>(new g4va(G4Colour(0.93,0.51,0.93))); //Violet
  LightYellowVisAtt = unique_ptr<g4va>(new g4va(G4Colour(0.8,0.8,0.2))); //kapton light yellow
  
  BlackVisAtt = unique_ptr<g4va>(new g4va(G4Colour(0.0,0.5,1.0)));	//black(dark blue)
  CableVisAtt = unique_ptr<g4va>(new g4va(G4Colour(1.0,1.0,1.0))); //white

  RedVisAtt = unique_ptr<g4va>(new g4va(G4Colour(1.0,0.0,0.0)));	//red
  YellowGreenVisAtt=unique_ptr<g4va>(new g4va(G4Colour(153./255.,204./255.,50./255.))); 
  LightPurpleVisAtt=unique_ptr<g4va>(new g4va(G4Colour(155./255.,48./255.,255./255.))); 
  PurpleVisAtt=unique_ptr<g4va>(new g4va(G4Colour(128./255.,0./255.,128./255.))); 
  DarkOrangeVisAtt=unique_ptr<g4va>(new g4va(G4Colour(255./255.,140./255.,0./255.))); 
  SkyBlueVisAtt=unique_ptr<g4va>(new g4va(G4Colour(0./255.,127./255.,255./255.))); 
  LightGreenVisAtt=unique_ptr<g4va>(new g4va(G4Colour(102/255.,153./255.,51./255.))); 
  LightBlueVisAtt=unique_ptr<g4va>(new g4va(G4Colour(132/255.,112./255.,250./255.))); 
  
  SteelVisAtt = unique_ptr<g4va>(new g4va(G4Colour(176./255.,196./255.,222./255.)));
  IronVisAtt = unique_ptr<g4va>(new g4va(G4Colour(100./255.,149./255.,237./255.)));
  SilverVisAtt = unique_ptr<g4va>(new g4va(G4Colour(153/255.,204./255.,255./255.)));  
  LeadVisAtt = unique_ptr<g4va>(new g4va(G4Colour(204/255.,204./255.,255./255.)));

  MagFieldVisAtt = unique_ptr<g4va>(new g4va(G4Colour(1., 0., 1.))); 
  
}

