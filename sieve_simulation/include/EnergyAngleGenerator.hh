#ifndef EnergyAngleGenerator_h 
#define EnergyAngleGenerator_h

#include <random>
#include <vector>

#include "Bound.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

namespace B1
{

// an implementation of a MHMC aglorithm using a table of energy-angle values
class EnergyAngleGenerator {
private: 

    G4ParticleDefinition* fParticleDef; 

    //table
    std::vector<double> fTable; 
    int fNpts_energy, fNpts_cosTheta; 

    //bounds for generation
    Bound<double> fBound_cosTheta, fBound_energy; 

    //current state of the generator
    double fCosTheta, fEnergy, fAmplitude; 

    double fScanRate_energy, fScanRate_cosTheta; 


    inline double Table(int iE, int iC) const { return fTable[iE*fNpts_cosTheta + iC]; }

    inline double RndmRange(double min, double max) const { return min + (max-min)*G4UniformRand(); }

public: 

    EnergyAngleGenerator(
        const std::vector<double>& _data,
        int npts_energy, int npts_cos_theta, 
        double min_energy, double max_energy, 
        double min_costheta, double max_costheta, 
        G4ParticleDefinition* particle
    ); 
    
    EnergyAngleGenerator() {};

    //return the amplitude at a particular phase-space point
    double Amplitude(double energy, double cos_theta) const; 

    //get particle definition
    G4ParticleDefinition* GetDefinition() { return fParticleDef; } 

    double GetEnergy() const { return fEnergy; }
    double GetCosTheta() const { return fCosTheta; }

    void SetEnergyScanRate(double _x) { fScanRate_energy=_x; } 
    void SetCosThetaScanRate(double _x) { fScanRate_cosTheta=_x; } 

    //do a metrpolis update 
    void Update(); 

};

}; 

#endif