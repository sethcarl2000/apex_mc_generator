#ifndef B1DetectorMessenger_h
#define B1DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4String.hh"
#include "G4UIcommand.hh"
//stdlib 
#include <memory> 
#include <functional> 
#include <vector> 


namespace B1 
{

class DetectorConstruction; 

class DetectorMessenger : public G4UImessenger {
public: 

    DetectorMessenger(DetectorConstruction*); 
    ~DetectorMessenger(); 

    void SetNewValue(G4UIcommand*, G4String); 

private: 

    DetectorConstruction* fDetector; 

    /// @brief basic struct to store a UI command, and a std::funciton which propagates this command to the DetectorConstruction Class 
    struct UIcommand_t {
        /// pointer to the UI command
        G4UIcommand* ptr; 

        /// funciton which propagates the information from this command to the detector
        std::function<void(G4UIcommand*,G4String)> fcn; 
    }; 

    /// A list of all UI commands.
    std::vector<UIcommand_t> fCommands; 

    /// @brief a small macro add a G4String-command to the list of commands 
    /// @param cmd_name name of the command, as it appears in a macro
    /// @param parameter_name name of the command parameter (i don't know what this is)
    /// @param signal_slot signal slot in target class to propagate this command to (see implementation in this class's constructor for an example)
    /// @param default_value default value of the command
    /// @param candidates list of valid command options, which should be given in a single string, with each option space-separated
    /// @param guidance command guidance
    void AddUIcommand_string(
        G4String cmd_name, 
        G4String parameter_name, 
        void (DetectorConstruction::*signal_slot)(G4String), 
        G4String default_value,
        G4String cadidates,
        G4String guidance=""
    ); 
}; 

}

#endif 