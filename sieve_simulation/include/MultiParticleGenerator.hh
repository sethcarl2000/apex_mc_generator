#ifndef MultiParticleGenerator_h
#define MultiParticleGenerator_h 1

#include <memory> 
#include <vector> 
#include "EnergyAngleGenerator.hh"
#include "G4ParticleDefinition.hh"

namespace B1 
{

class MultiParticleGenerator { 
private: 
    //'stack' of probabilities (must add to 1)
    std::vector<double> fProbStack;

    std::vector<EnergyAngleGenerator*> fGenerators; 

    //current state of generators
    double fCosTheta, fEnergy; 
    G4ParticleDefinition *fParticleDef; 

    //update 'probability stack' 
    void UpdateProbStack(); 

    G4int fVerbose; 

public: 

    MultiParticleGenerator();
    ~MultiParticleGenerator();

    /// @brief update the generator
    void Update(); 

    /// @brief add process to the list of processes  
    void AddProcess(EnergyAngleGenerator* process); 

    double GetCosTheta() const { return fCosTheta; } 

    double GetEnergy()   const { return fEnergy; }

    G4ParticleDefinition* GetParticleDef() const { return fParticleDef; }
};

};

#endif