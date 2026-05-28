

#include "MultiParticleGenerator.hh"
#include "Randomize.hh"

#include "RunParameters.hh"

#include <algorithm> 
#include <cstdio> 

namespace B1
{
//____________________________________________________________________________________________
MultiParticleGenerator::~MultiParticleGenerator()
{
    for (auto& process : fGenerators) { if (process) delete process; process = nullptr; } 
    fGenerators.clear(); 
}
//____________________________________________________________________________________________

//____________________________________________________________________________________________
void MultiParticleGenerator::AddProcess(EnergyAngleGenerator* process)
{
    //sort all cross sections so that the largest cross-sections are first in the list. 
    for (auto it = fGenerators.begin(); it != fGenerators.end(); it++) {

        if (process->GetTotalCS() > (*it)->GetTotalCS()) {

            fGenerators.insert(it, process); 
            break; 
        }
    }

    //update the list of process probabilities 
    UpdateProbStack(); 

    if (RunParameters::Instance()->Verbosity_generator() >= 3) {
        std::printf(
            " in <MultiParticleGenerator>:  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
            "   new process added. total: %zi\n"
            "   processes: \n",
            fGenerators.size()
        );

        double last_prob=0.; 
        for (size_t i=0; i<fGenerators.size(); i++) {

            const auto process = fGenerators[i];
            std::printf(
                "   --- %zi/%zi, cross-section: % .5e GeV^-2, (%5.1f%%). process: %s\n",
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
    //pick a random particle to generate
    double p = G4UniformRand(); 
    int i=0; 
    while (p < fProbStack.at(i)) { 
        std::printf("p[%i] = % .8e\n",i,fProbStack.at(i)); 
        ++i; 
    }
    
    std::cout << "Index is " << i << std::endl; 

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