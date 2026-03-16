
#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4ApplicationState.hh"

//specific UI commands 
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"

using namespace std; 


namespace B1 
{

//_____________________________________________________________________________________________________________
DetectorMessenger::DetectorMessenger(DetectorConstruction* detector)
    : fDetector{detector}
{
    G4String command_prefix = "/detector/"; 

    /// add all the commands we want here ---------------------------------------------------------
    
    // Command to select which arm to use
    AddUIcommand_string(
        command_prefix + "arm",  //command name
        "HRS_arm",               // parameter name
        &DetectorConstruction::SetArm, // associated function in DetectorConstruction
        "RHRS",                  // default value
        "RHRS LHRS",             // possible valid inputs
        "Deterimes whether the LHRS or RHRS will be simulated. valid inputs are 'RHRS or LHRS'"
    );
}
//_____________________________________________________________________________________________________________
DetectorMessenger::~DetectorMessenger()
{
    /// delete all commands
    for (auto& command : fCommands) {
        if (command.ptr != nullptr) delete command.ptr;  
    }
}
//_____________________________________________________________________________________________________________
void DetectorMessenger::AddUIcommand_string(
    G4String command_name, 
    G4String parameter_name, 
    void (DetectorConstruction::*signal_slot)(G4String), 
    G4String default_value, 
    G4String candidates, 
    G4String guidance
)
{ 
    auto cmd = new G4UIcmdWithAString(command_name, this);

    cmd->SetParameterName(parameter_name, false); 
    cmd->SetDefaultValue(default_value);
    cmd->SetCandidates(candidates); 
    cmd->SetGuidance(guidance); 

    /// construct a 'signal function' which tells the DetectorConstruction class about our updated command
    auto signal_fcn = [this, signal_slot](G4UIcommand* parent_ptr, G4String new_val)
    {
        auto cmd_ptr = dynamic_cast<G4UIcmdWithAString*>(parent_ptr);
        (fDetector->*signal_slot)(new_val);
        return; 
    }; 
    
    // our new UI command  
    UIcommand_t new_ui_command {
        .ptr = cmd, 
        .fcn = signal_fcn
    }; 
    
    //add this command to the list of all commands
    fCommands.push_back(new_ui_command); 
}
//_____________________________________________________________________________________________________________
void DetectorMessenger::SetNewValue(G4UIcommand* issued_command, G4String value)
{
    //scan our list of commands to see which was issued, then execute the issued command
    for (auto& ui_command : fCommands) {

        if (issued_command == ui_command.ptr) { ui_command.fcn(issued_command, value); } 
    }
}
//_____________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________

}