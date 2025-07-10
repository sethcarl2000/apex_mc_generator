// ********************************************************************
// $Id: HRSVisAttribute.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//

#ifndef HRSVisAttribute_h
#define HRSVisAttribute_h 1

#include <memory> 

class G4VisAttributes;

class HRSVisAttribute
{
public:
  HRSVisAttribute();
  virtual ~HRSVisAttribute();

private:
  void ConstructVisAttribute();

public:
  
  std::unique_ptr<G4VisAttributes> Invisible_VisAtt; 
  
  std::unique_ptr<G4VisAttributes> HallVisAtt;
  std::unique_ptr<G4VisAttributes> MagneticVisAtt;
  std::unique_ptr<G4VisAttributes> WhiteVisAtt;
  std::unique_ptr<G4VisAttributes> DarkRedVisAtt;
  std::unique_ptr<G4VisAttributes> OrangeVisAtt;
  std::unique_ptr<G4VisAttributes> GrayVisAtt;
  std::unique_ptr<G4VisAttributes> YellowVisAtt;
  std::unique_ptr<G4VisAttributes> PcbGreenVisAtt;
  std::unique_ptr<G4VisAttributes> CuBrownVisAtt;
  std::unique_ptr<G4VisAttributes> DarkBlueVisAtt;		//36;24;130
  std::unique_ptr<G4VisAttributes> VioletVisAtt;			//238;130;238
  std::unique_ptr<G4VisAttributes> LightYellowVisAtt;

  std::unique_ptr<G4VisAttributes> BlackVisAtt; 
  std::unique_ptr<G4VisAttributes> CableVisAtt;

  std::unique_ptr<G4VisAttributes> RedVisAtt;				//255;204;50
  std::unique_ptr<G4VisAttributes> YellowGreenVisAtt;		//153;204;50
  std::unique_ptr<G4VisAttributes> LightPurpleVisAtt;		//155;48;255
  std::unique_ptr<G4VisAttributes> PurpleVisAtt;			//128,0,128
  std::unique_ptr<G4VisAttributes> DarkOrangeVisAtt;		//255,140,0
  std::unique_ptr<G4VisAttributes> SkyBlueVisAtt;			//0;127;255
  std::unique_ptr<G4VisAttributes> LightGreenVisAtt;
  std::unique_ptr<G4VisAttributes> LightBlueVisAtt;

  std::unique_ptr<G4VisAttributes> SteelVisAtt;
  std::unique_ptr<G4VisAttributes> IronVisAtt;
  std::unique_ptr<G4VisAttributes> SilverVisAtt;
  std::unique_ptr<G4VisAttributes> LeadVisAtt;

  std::unique_ptr<G4VisAttributes> MagFieldVisAtt; 
  
};

#endif

