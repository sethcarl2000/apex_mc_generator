

#include "MultiParticleGenerator.hh"
#include "Randomize.hh"

#include "G4Exception.hh"

#include "RunParameters.hh"

#include <algorithm> 
#include <cstdio> 

namespace B1
{
//____________________________________________________________________________________________
//____________________________________________________________________________________________
//____________________________________________________________________________________________
MultiParticleGenerator::MultiParticleGenerator()
    : fGenerators{}, fProbStack{}
{
    fVerbose = RunParameters::Instance()->Verbosity_generator(); 
}
//____________________________________________________________________________________________
MultiParticleGenerator::~MultiParticleGenerator()
{
    for (auto& process : fGenerators) { if (process) delete process; process = nullptr; } 
    fGenerators.clear(); 
}
//____________________________________________________________________________________________
void MultiParticleGenerator::AddProcess(EnergyAngleGenerator* process)
{
    //sort all cross sections so that the largest cross-sections are first in the list. 
    size_t i=0; 
    while (i < fGenerators.size() && (process->GetTotalCS() < fGenerators[i]->GetTotalCS())) { ++i; }

    fGenerators.insert( fGenerators.begin()+i, process );

    //update the list of process probabilities 
    UpdateProbStack(); 

    if (fVerbose >= 3) {
        std::printf(
            " in <MultiParticleGenerator>:  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
            "   new process added: '%s'\n"
            "   total n. processes: %zi\n"
            "   processes: \n",
            process->GetDescription().c_str(),
            fGenerators.size()
        );

        double last_prob=0.; 
        for (size_t i=0; i<fGenerators.size(); i++) {

            const auto process = fGenerators[i];
            std::printf(
                "   --- %zi/%zi, cross-section: % .5e MeV^-2, (%7.4f%%). process: %s\n",
                i,fGenerators.size(), 
                process->GetTotalCS(), 
                100.*(fProbStack[i]-last_prob), 
                process->GetDescription().c_str()
            );
            last_prob = fProbStack[i];
        }
        
        std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    }

}
//____________________________________________________________________________________________
void MultiParticleGenerator::UpdateProbStack()
{
    fProbStack.clear(); 
    fProbStack.reserve(fGenerators.size());

    double total_cs=0.; 

    for (size_t i=0; i<fGenerators.size(); i++) {

        const auto& process = fGenerators[i]; 

        total_cs = total_cs + process->GetTotalCS(); 

        fProbStack.push_back( total_cs );
    }; 

    for (auto& p : fProbStack) { p = p/total_cs; }
}
//____________________________________________________________________________________________
void MultiParticleGenerator::Update()
{
    if (fGenerators.empty()) {
        G4Exception(
            "MultiParticleGenerator::Update()",
            "No processes have been added before Update() was called", 
            G4ExceptionSeverity::RunMustBeAborted, 
            "'Update()' cannot be called when no processes have been added."
        );
        return; 
    }

    //pick a random particle to generate
    double p = G4UniformRand(); 
    int i=0; 
    while (p > fProbStack[i]) { 
        if (fVerbose >= 3) { std::printf("P[%i] = % .8f\n", i, fProbStack.at(i)); } 
        ++i; 
    }
    
    if (i >= fGenerators.size()) {
        G4Exception(
            "MutliParticleGenerator::Update()", 
            "Invalid process index",
            G4ExceptionSeverity::RunMustBeAborted,
            "An invalid index was reached trying to generate a process."
        ); 
        return; 
    }
    auto& process = fGenerators[i];

    process->Update(); 

    fParticleDef = process->GetDefinition(); 
    fCosTheta = process->GetCosTheta();
    fEnergy   = process->GetEnergy(); 
}
//____________________________________________________________________________________________
//____________________________________________________________________________________________
//____________________________________________________________________________________________
//____________________________________________________________________________________________
//____________________________________________________________________________________________
//____________________________________________________________________________________________

};