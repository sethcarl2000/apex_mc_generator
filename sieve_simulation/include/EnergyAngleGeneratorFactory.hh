#ifndef EnergyAngleGeneratorFactory_h
#define EnergyAngleGeneratorFactory_h 

#include "EnergyAngleGenerator.hh"
#include "G4String.hh"
//this class will instantiate an energy-angle generator, when requested. 

namespace B1 
{

class EnergyAngleGeneratorFactory {
private: 

    EnergyAngleGeneratorFactory() {};
    ~EnergyAngleGeneratorFactory() {}; 

public: 
    
    /// @return Energy-Angle generator corresponding to electron elastic scattering
    static EnergyAngleGenerator* Elastic(); 

    /// @return Energy-Angle generator corresponding to electron scattering off nuclear target (with photon produced)
    static EnergyAngleGenerator* BetheHeitlerPhotoproduction(); 

    /// @return Energy-Angle generator corresponding to electron production in trident process
    static EnergyAngleGenerator* Trident_Electron(); 

    /// @return Energy-Angle generator corresponding to positron production in trident process
    static EnergyAngleGenerator* Trident_Positron(); 
};

};

#endif