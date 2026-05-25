#include <EnergyAngleGenerator.hh> 
#include <cmath> 
#include <stdio.h>
#include "Randomize.hh"

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

    fScanRate_cosTheta = fBound_cosTheta.span()/2.;
    fScanRate_energy   = 25.; 
}; 
//_____________________________________________________________________________________________
double EnergyAngleGenerator::Amplitude(double energy, double cos_theta) const
{
    double dE = fBound_energy.span()/((double)fNpts_energy-1);
    double dC = fBound_cosTheta.span()/((double)fNpts_cosTheta-1); 

    Bound<int> bound_ind_E{ 0, fNpts_energy-1 };
    Bound<int> bound_ind_C{ 0, fNpts_cosTheta-1 };

    int iE_lo = std::floor( (energy - fBound_energy.min)/dE ); 
    int iC_lo = std::floor( (cos_theta - fBound_cosTheta.min)/dC ); 

    double iE_hi = iE_lo+1; 
    double iC_hi = iC_lo+1; 

    iE_lo = bound_ind_E.enforce(iE_lo); 
    iE_hi = bound_ind_E.enforce(iE_hi); 
    
    iC_lo = bound_ind_C.enforce(iC_lo); 
    iC_hi = bound_ind_C.enforce(iC_hi); 

    //first, interpolate with the energy 
    double E0 = fBound_energy.min + dE*((double)iE_lo); 
    double E1 = E0 + dE; 

    double C0 = fBound_cosTheta.min + dC*((double)iC_lo);
    double C1 = C0 + dC; 

    double delta_E = (energy - E0)/dE; 
    double delta_C = (cos_theta - C0)/dC; 

    return
        (1-delta_C)*( Table(iE_lo,iC_lo) + delta_E*(Table(iE_hi,iC_lo) - Table(iE_lo,iC_lo)) ); 
        + (delta_C)*( Table(iE_lo,iC_hi) + delta_E*(Table(iE_hi,iC_hi) - Table(iE_lo,iC_hi)) );     
}
//_____________________________________________________________________________________________
void EnergyAngleGenerator::Update()
{
    //propose update 
    //fCosTheta += (1. - 2.*Rndm())*fScanRate_cosTheta; 

    double cos_theta = fBound_cosTheta.enforce( 
        fCosTheta + (1. - 2.*G4UniformRand())*fScanRate_cosTheta 
    );

    double energy = fBound_energy.enforce( 
        fEnergy + (1. - 2.*G4UniformRand())*fScanRate_energy
    );

    double new_amplitude = Amplitude(energy, cos_theta);

    std::printf("energy: %6.1f MeV, cos(theta): %6.4f, new / old amplitude: %.3f;",
        energy, cos_theta, new_amplitude/fAmplitude
    );

    if ( G4UniformRand() < new_amplitude/fAmplitude ) {
        std::printf(" accepted\n");
        //accept update 
        fEnergy = energy;
        fCosTheta = cos_theta; 
        fAmplitude = new_amplitude;
    } else {
        std::printf("\n"); 
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