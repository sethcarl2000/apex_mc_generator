#ifndef UserMessenger_h
#define UserMessenger_h 1

#include "G4UIcommand.hh"
#include "G4UImessenger.hh"
#include "G4String.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include <functional>
#include <vector>

namespace B1 
{

/// @brief This class can generate different G4UIcommands for the target class 
/// @tparam UserClass 
template<typename UserClass> class UserMessenger : public G4UImessenger {
private: 
    
    /// @brief basic struct to store a UI command, and a std::funciton which propagates this command to the DetectorConstruction Class 
    struct UIcommand_t {
        /// pointer to the UI command
        G4UIcommand* ptr; 

        /// funciton which propagates the information from this command to the detector
        std::function<void(G4UIcommand*,G4String)> fcn; 
    }; 

    //the instance of the target class to point commands at
    UserClass *fTarget;  

    std::vector<UIcommand_t> fCommands; 

public: 
    UserMessenger(UserClass *target) : fTarget{target}, fCommands{} {}; 

    ~UserMessenger() { 
        for (auto& command : fCommands) if (command.ptr) delete command.ptr; 
    };

    /// @brief Set new value of a parameter using a UI command  
    void SetNewValue(G4UIcommand* issued_command, G4String new_val) {

        for (auto& UI_command : fCommands) {
            if (UI_command.ptr == issued_command) { UI_command.fcn(issued_command,new_val); }
        }
    }; 

    /// @brief a small macro add a G4String-command to the list of commands 
    /// @param cmd_name name of the command, as it appears in a macro
    /// @param parameter_name name of the command parameter (i don't know what this is)
    /// @param signal_slot signal slot in target class to propagate this command to (see implementation in this class's constructor for an example)
    /// @param default_value default value of the command
    /// @param candidates list of valid command options, which should be given in a single string, with each option space-separated
    /// @param guidance command guidance
    void AddCommand_string(
        G4String cmd_name, 
        G4String parameter_name, 
        void (UserClass::*signal_slot)(G4String), 
        G4String default_value,
        G4String candidates,
        G4String guidance=""
    )
    { 
        auto cmd = new G4UIcmdWithAString(cmd_name, this);

        cmd->SetParameterName(parameter_name, this); 
        cmd->SetDefaultValue(default_value);
        if (candidates != "") cmd->SetCandidates(candidates); 
        if (guidance != "")   cmd->SetGuidance(guidance); 

        /// construct a 'signal function' which tells the DetectorConstruction class about our updated command
        auto signal_fcn = [this, signal_slot](G4UIcommand* parent_ptr, G4String new_val)
        {
            auto cmd_ptr = dynamic_cast<G4UIcmdWithAString*>(parent_ptr);
            (fTarget->*signal_slot)(new_val);
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

    /// @brief a small macro add a double-with-unit command to the list of commands 
    /// @param cmd_name name of the command, as it appears in a macro
    /// @param parameter_name name of the command parameter (i don't know what this is)
    /// @param signal_slot signal slot in target class to propagate this command to (see implementation in this class's constructor for an example)
    /// @param default_value default value of the command
    /// @param default_unit default unit of quantity 
    /// @param guidance command guidance
    void AddCommand_doubleWithUnit(
        G4String cmd_name, 
        G4String parameter_name, 
        void (UserClass::*signal_slot)(G4double), 
        G4double default_value,
        G4String default_unit,
        G4String guidance=""
    )
    { 
        auto cmd = new G4UIcmdWithADoubleAndUnit(cmd_name, this);

        cmd->SetParameterName(parameter_name, this); 
        cmd->SetDefaultValue(default_value);
        cmd->SetDefaultUnit(default_unit); 
        if (guidance != "")   cmd->SetGuidance(guidance); 

        /// construct a 'signal function' which tells the DetectorConstruction class about our updated command
        auto signal_fcn = [this, signal_slot](G4UIcommand* parent_ptr, G4String new_val)
        {
            auto cmd_ptr = dynamic_cast<G4UIcmdWithADoubleAndUnit*>(parent_ptr);
            (fTarget->*signal_slot)(cmd_ptr->GetNewDoubleValue(new_val));
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

};

}

#endif 