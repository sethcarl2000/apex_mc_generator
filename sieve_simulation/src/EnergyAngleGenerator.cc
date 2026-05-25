#include <EnergyAngleGenerator.hh> 
#include <cmath> 
#include <stdio.h>
#include <limits> 
#include "Randomize.hh"

namespace {
    constexpr double kInfinity = std::numeric_limits<double>::infinity(); 
};

namespace B1
{

//_____________________________________________________________________________________________
EnergyAngleGenerator::EnergyAngleGenerator(
    const std::vector<double>& _data,
    int npts_energy, int npts_cos_theta, 
    double min_energy, double max_energy, 
    double min_costheta, double max_costheta, 
    G4ParticleDefinition* particle
)
    : fTable{_data}, 
    fNpts_energy{npts_energy},
    fNpts_cosTheta{npts_cos_theta},
    fBound_energy{min_energy, max_energy}, 
    fBound_cosTheta{min_costheta, max_costheta}, 
    fParticleDef{particle}
{
    //start with the largest value of Cos & Energy 
    fCosTheta   = fBound_cosTheta.max; 
    fEnergy     = fBound_energy.max; 
    fAmplitude  = Amplitude(fEnergy, fCosTheta);

    fScanRate_cosTheta = fBound_cosTheta.span()/8.;
    fScanRate_energy   = fBound_energy.span()/100.; 
}; 
//_____________________________________________________________________________________________
double EnergyAngleGenerator::Amplitude(double energy, double cos_theta) const
{
    bool single_energy_pt = fNpts_energy==1;
    bool single_cos_pt    = fNpts_cosTheta==1;   

    double dE = single_energy_pt ? kInfinity : fBound_energy.span()/((double)fNpts_energy-1);
    double dC = single_cos_pt    ? kInfinity : fBound_cosTheta.span()/((double)fNpts_cosTheta-1); 

    Bound<int> bound_ind_E{ 0, fNpts_energy-1 };
    Bound<int> bound_ind_C{ 0, fNpts_cosTheta-1 };

    int iE_lo = std::floor( (energy - fBound_energy.min)/dE ); 
    int iC_lo = std::floor( (cos_theta - fBound_cosTheta.min)/dC ); 

    int iE_hi, iC_hi; 

    iE_lo = bound_ind_E.enforce(iE_lo); 
    iE_hi = bound_ind_E.enforce(iE_lo+1); 
    
    iC_lo = bound_ind_C.enforce(iC_lo); 
    iC_hi = bound_ind_C.enforce(iC_lo+1); 

    //first, interpolate with the energy 
    double E0, E1;
    
    if (single_energy_pt) {
        E0 = E1 = fBound_energy.min;
    } else {
        E0 = fBound_energy.min + dE*((double)iE_lo); 
        E1 = E0 + dE; 
    }
    
    double C0, C1; 

    if (single_cos_pt) {
        C0 = C1 = fBound_cosTheta.min; 
    } else {
        C0 = fBound_cosTheta.min + dC*((double)iC_lo);
        C1 = C0 + dC; 
    }

    double delta_E = single_energy_pt ? 0. : (energy - E0)/dE; 
    double delta_C = single_cos_pt ? 0. : (cos_theta - C0)/dC; 

    double Amp_lo = Table(iE_lo,iC_lo) + delta_E*(Table(iE_hi,iC_lo) - Table(iE_lo,iC_lo));
    double Amp_hi = Table(iE_lo,iC_hi) + delta_E*(Table(iE_hi,iC_hi) - Table(iE_lo,iC_hi));

    return (1.-delta_C)*( Amp_lo ) + (delta_C)*( Amp_hi );     
}
//_____________________________________________________________________________________________
void EnergyAngleGenerator::Update()
{
    //propose update 
    //fCosTheta += (1. - 2.*Rndm())*fScanRate_cosTheta; 

    double cos_theta_min = fBound_cosTheta.enforce( fCosTheta - fScanRate_cosTheta );
    double cos_theta_max = fBound_cosTheta.enforce( fCosTheta + fScanRate_cosTheta );

    double energy_min = fBound_energy.enforce( fEnergy - fScanRate_energy );
    double energy_max = fBound_energy.enforce( fEnergy + fScanRate_energy );

    double energy = RndmRange( energy_min, energy_max );
    double cos_theta = RndmRange( cos_theta_min, cos_theta_max ); 

    double new_amplitude = Amplitude(energy, cos_theta);

    /*std::printf("energy: %6.1f MeV, cos(theta): %6.4f, new / old amplitude: %.3f;",
        energy, cos_theta, new_amplitude/fAmplitude
    );*/

    if ( G4UniformRand() < new_amplitude/fAmplitude ) {
        //std::printf(" accepted\n");
        //accept update 
        fEnergy = energy;
        fCosTheta = cos_theta; 
        fAmplitude = new_amplitude;
    } else {
        //std::printf("\n"); 
    }
}
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________

};