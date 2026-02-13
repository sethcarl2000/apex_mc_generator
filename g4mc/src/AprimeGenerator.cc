
#include "AprimeGenerator.hh"
#include "ApexTargetGeometry.hh"
#include <stdio.h>

using namespace std;

namespace {

  constexpr double kMass_e = 0.501; ///Electron mass, MeV/c^2
  
  constexpr int kN_foils = 10;      /// number of production foils  
  
  /// ---- parameters for search of phase-space neighborhood
  constexpr double kX_sweep_range     = 2.0e-5;    //maximum change in x
  constexpr double kTheta_sweep_range = 0.300e-3;  //maximum change in A' theta
  constexpr double kPhi_sweep_range   = CLHEP::pi/18.;    //maximum change in A' phi
  constexpr double kFoil_jump_prob    = 0.20;      //probablility to 'jump' from the current foil\
 to another foil
  constexpr double kRotation_RMS      = 0.040/1.73205;  /// 0.020 is the average magnitude of the random rotation we perform on the electron's A'-rest-frame orientation, in radians (in A' rest frame).

}


//__________________________________________________________________________________________
AprimeGenerator::AprimeGenerator(double _mass_A, double _beam_E, double _spec_p0)
    : fMass_Aprime{_mass_A},
    fBeam_energy{_beam_E},
    fSpectrometer_p0{_spec_p0},
    fE_electron{fMass_Aprime/2.}
{
    std::random_device rd;

    fTwister = std::mt19937(rd()); 

#ifdef DEBUG
    printf("Aprime constructor ~~~~~~~ \n"
	   " m_A      %.f MeV/c^2\n"
	   " E_e+/-   %.f MeV\n",
	   fMass_Aprime, fE_electron);
#endif
    
    ResetState(); 
}
//__________________________________________________________________________________________
void AprimeGenerator::ResetState()
{
  fE_electron = fMass_Aprime/2.;
  
  fState = State_t{
    .vertex     = G4ThreeVector(0., 0., Get_foil_z(4)), 
    .x          = 1.,
    .theta      = 0.,
    .phi        = 0., 
    .P_electron = G4ThreeVector(sqrt(fE_electron*fE_electron - kMass_e*kMass_e), 0., 0.),
    .foil_num   = 4, 
    .amplitude  = A_production_amplitude(1., 0.), /// a negative amplitude indicates that this state is invalid  
    .in_acceptance=false
  };
}
//__________________________________________________________________________________________
G4ThreeVector AprimeGenerator::Random_rotation(const G4ThreeVector& v, double RMS_rotation_mag) 
{
    //implementing: vR = v + v x R 
    //              vR_i = \delta_ij v_j   +  \epsilon_ijk R_j v_k    
    // where: 
    //      \delta_ij is the kroneker delta tensor, and 
    //      \epsilon_ijk is the antisymmetric tensor. 
    // and we use the fact that 'R' should be a small vector, so this is a first-order approximation w/r/t mag(R). 
    double vx(v[0]), vy(v[1]), vz(v[2]); 
    G4ThreeVector vR( vx, vy, vz );  

    double mag = v.mag(); 

    double RX = RMS_rotation_mag*Gaus();     
    double RY = RMS_rotation_mag*Gaus();
    double RZ = RMS_rotation_mag*Gaus();
    
    vR[0] +=  RY*vz - RZ*vy;  
    vR[1] +=  RZ*vx - RX*vz; 
    vR[2] +=  RX*vy - RY*vx; 

    return vR.unit() * mag; 
}   
//__________________________________________________________________________________________
double AprimeGenerator::A_production_amplitude(double x, double theta_A)
{
    if (x < 0. || x > 1.) return 0.;

    double th_A2   = theta_A*theta_A; 
    double beam_E2 = fBeam_energy*fBeam_energy; 
    double m_A2    = fMass_Aprime*fMass_Aprime; 

    double U = beam_E2 * x * th_A2   +   m_A2*(1. - x)/x   +  (kMass_e*kMass_e)*x; 

    return beam_E2*x*(  1. - x + (x*x/2.) - x*x*(1-x)*m_A2*beam_E2*th_A2/(U*U)  )/(U*U); 
}
//__________________________________________________________________________________________
void AprimeGenerator::SetRange_x_sv(bool is_RHRS, double min, double max)
{
    if (is_RHRS) { fRange_R_x_sv[0] =min; fRange_R_x_sv[1] =max; }
    else         { fRange_L_x_sv[0] =min; fRange_L_x_sv[1] =max; }
}
//__________________________________________________________________________________________
void AprimeGenerator::SetRange_y_sv(bool is_RHRS, double min, double max)
{
    if (is_RHRS) { fRange_R_y_sv[0] =min; fRange_R_y_sv[1] =max; }
    else         { fRange_L_y_sv[0] =min; fRange_L_y_sv[1] =max; }
}
//__________________________________________________________________________________________
void AprimeGenerator::Update()
{
  fN_updates++; 
  
    auto new_point = fState; 

    new_point.x     += (1. - 2.*Rand())*kX_sweep_range; 
    new_point.theta += (1. - 2.*Rand())*kTheta_sweep_range; 
    new_point.phi   += (1. - 2.*Rand())*kPhi_sweep_range; 

    double jump_r = Rand(); 
    if (jump_r < kFoil_jump_prob) {
        //decide whether to jump to the next or previous foil
        if (jump_r / kFoil_jump_prob < 0.5) { new_point.foil_num += +1; } else { new_point.foil_num += -1; } 
        //fix the new foil num if its out-of-range
        if (new_point.foil_num < 0) new_point.foil_num = 0; 
        if (new_point.foil_num >= kN_foils) new_point.foil_num = kN_foils-1; 
    }
    
    //get a new electorn direction (in A' rest frame)
    new_point.P_electron = Random_rotation(new_point.P_electron, kRotation_RMS);  

    double p_electron_mag = std::sqrt( fE_electron*fE_electron - kMass_e*kMass_e ); 

    new_point.P_electron = new_point.P_electron.unit() * p_electron_mag; 

    double E_a = new_point.x * fBeam_energy; 
    double p_A_mag = sqrt( (E_a*E_a) - fMass_Aprime*fMass_Aprime ); 

    G4ThreeVector p_A(
        p_A_mag * cos(new_point.theta) * cos(new_point.phi), 
        p_A_mag * sin(new_point.theta) * sin(new_point.phi), 
        p_A_mag * cos(new_point.theta)
    ); 

    //now, in the rest-frame of the A', we generate the directions of the positron / electron randomly
    double gamma_A = E_a / fMass_Aprime; 
    double beta_A  = sqrt( 1. - 1./(gamma_A*gamma_A) ); 

    //______________________________________________________________________________________
    //take a vector in the rest-frame of the A', and boost it to the lab frame
    auto boost_and_rotate = [gamma_A, beta_A, &new_point](const G4ThreeVector &p, double E)
    {
        G4ThreeVector p_boost(
            p[0], 
            p[1], 
            gamma_A * (p[2] + beta_A*E)
        ); 
        
        //so, now that we've boosed the momentum back to the lab, 
        //let's rotate this vector, so that it's boost is in line with the momentum of the A' 
        p_boost.rotateX( new_point.theta );
        p_boost.rotateZ( new_point.phi );
        
        return p_boost; 
    };
    //_____________________________________________________________________________________
    
    //randomly generate the vertex, in line with the production target geometry        
    new_point.vertex = G4ThreeVector( 
        2.*( 1. - 2.*Rand() ), 
        2.*( 1. - 2.*Rand() ), 
        Get_foil_z(new_point.foil_num)
    );
    
    //check the acceptances of the electron and positron
    new_point.Pe_boosted = boost_and_rotate(new_point.P_electron, fE_electron); 
    new_point.Pp_boosted = boost_and_rotate(-1.*new_point.P_electron, fE_electron); 

    new_point.in_acceptance = true; 
    
    //reject update if the new electron is outside the momentum range
    if (fabs(new_point.Pe_boosted.mag() - fSpectrometer_p0)/fSpectrometer_p0 > fRange_dp) { new_point.in_acceptance=false; } 
    if (fabs(new_point.Pp_boosted.mag() - fSpectrometer_p0)/fSpectrometer_p0 > fRange_dp) { new_point.in_acceptance=false; } 

#ifdef DEBUG
    
    double L_x_sv, L_y_sv, R_x_sv, R_y_sv; 

    Project_onto_sieve(true,  new_point.vertex, new_point.Pp_boosted, R_x_sv,R_y_sv);  
    Project_onto_sieve(false, new_point.vertex, new_point.Pe_boosted, L_x_sv,L_y_sv);  
        
    printf(
        " ~~~~~~~~~~~~~~~~ phase-space: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
	"   E_beam      %-.1f MeV\n"
	"   spec. p0    %-.1f MeV/c\n"
        "   x_hcs       %-3.4f mm\n"
        "   y_hcs       %-3.4f mm\n"
        "   foil #      %i\n"
        "   P[electron] - A' rest frame:  (%-.1f, %-.1f, %-.1f) MeV/c\n"
	"   A':\n"
	"       EA      = %-.1f\n"
	"       gamma   = %-.1f\n"
        "       EA/E0   = 1.0 - %.2e\n"
	"       theta_A = %-.3f mrad\n"
        "       phi_A   = %-3.1f deg\n"
        "   Electron:\n"
        "       x_sv = %-.3f mm    range [%-.3f, %-.3f]\n"
        "       y_sv = %-.3f mm    range [%-.3f, %-.3f]\n"
        "       P - lab: (%-.1f, %-.1f, %-.1f) MeV/c; E = %-.1f MeV.  range: %-.1f +/- 5\%\n"
        "   Positron:\n"
        "       x_sv = %-.3f mm    range [%-.3f, %-.3f]\n"
        "       y_sv = %-.3f mm    range [%-.3f, %-.3f]\n"
        "       P - lab: (%-.1f, %-.1f, %-.1f) MeV/c; E = %-.1f MeV.  range: %-.1f +/- 5\%\n"
        " inside acceptance? %s\n",
	fBeam_energy,
	fSpectrometer_p0,
        new_point.vertex.x(),
        new_point.vertex.y(), 
        new_point.foil_num,
        new_point.P_electron.x(), new_point.P_electron.y(), new_point.P_electron.z(),
	gamma_A, 
	E_a, 
	1. - new_point.x, 
        new_point.theta * 1.e3, 
        new_point.phi * (180./CLHEP::pi), 
        L_x_sv, fRange_L_x_sv[0], fRange_L_x_sv[1], 
        L_y_sv, fRange_L_y_sv[0], fRange_L_y_sv[1],
	new_point.Pe_boosted[0],new_point.Pe_boosted[1],new_point.Pe_boosted[2],
	new_point.Pe_boosted.mag(), fSpectrometer_p0, 
	R_x_sv, fRange_R_x_sv[0], fRange_R_x_sv[1],
        R_y_sv, fRange_R_y_sv[0], fRange_R_y_sv[1],
        new_point.Pp_boosted[0],new_point.Pp_boosted[1],new_point.Pp_boosted[2],
	new_point.Pp_boosted.mag(), fSpectrometer_p0,
	new_point.in_acceptance ? "yes" : "no"
    );
#endif 

    //if we're out of the acceptance w/r/t momentum, there's no reason to keep going.
    if (new_point.in_acceptance==false) return; 
    //_____________________________________________________________________________________
    auto check_acceptance = [&](const G4ThreeVector& p, bool is_RHRS)
    {   
        double x_sv, y_sv;  

	Project_onto_sieve(is_RHRS,
			   new_point.vertex,
			   p,
			   x_sv, y_sv);  
        
        double xmin = is_RHRS ? fRange_R_x_sv[0] : fRange_L_x_sv[0]; 
        double xmax = is_RHRS ? fRange_R_x_sv[1] : fRange_L_x_sv[1]; 

        if (x_sv > xmax || x_sv < xmin) return false;   
        
        double ymin = is_RHRS ? fRange_R_y_sv[0] : fRange_L_y_sv[0]; 
        double ymax = is_RHRS ? fRange_R_y_sv[1] : fRange_L_y_sv[1];     

        if (y_sv > ymax || y_sv < ymin) return false; 

        return true; 
    };
    //_____________________________________________________________________________________
    

    
    if (check_acceptance(new_point.Pp_boosted, true)  == false) { new_point.in_acceptance=false; } 
    if (check_acceptance(new_point.Pe_boosted, false) == false) { new_point.in_acceptance=false; } 


    //reject update
    if ((!new_point.in_acceptance) && (fState.in_acceptance)) { return; }
    
    new_point.amplitude = A_production_amplitude(new_point.x, new_point.theta);  

    //if we've gotten here, we can consider an update. 
    if ( new_point.amplitude > fState.amplitude || new_point.amplitude/fState.amplitude > Rand() ) {

        fState = new_point; 
    }
    return; 
}; 
//_________________________________________________________________________________________
void AprimeGenerator::Project_onto_sieve(bool is_RHRS, const G4ThreeVector& R1_hcs, const G4ThreeVector& dR_hcs, double &x_sv, double &y_sv)
{
  //this assumes that 'R' and 'dR' are both given in Hall Coordinate System (HCS)

  //R1 is a displacement vector, and dR is only a direction.
  //so we need to construct R2, a second displacement vector. 
  auto R2_hcs = R1_hcs + dR_hcs;

  auto R1 = ApexTargetGeometry::HCS_to_SCS(R1_hcs, is_RHRS); 
  auto R2 = ApexTargetGeometry::HCS_to_SCS(R2_hcs, is_RHRS);   
  
  auto S  = R2 - R1;

  R1 += -1.*ApexTargetGeometry::Get_sieve_pos(is_RHRS); 

  x_sv = R1.x() + (S.x()/S.z())*(0. - R1.z());
  y_sv = R1.y() + (S.y()/S.z())*(0. - R1.z());
  return; 
}
//_________________________________________________________________________________________
//_________________________________________________________________________________________
//_________________________________________________________________________________________
