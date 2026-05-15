#ifndef RunParameters_HH
#define RunParameters_HH

#include "UserMessenger.hh"
#include "ArmMode.hh"
#include "TargetMode.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"

/// @brief Meyers-Singleton class which stores global run parameters, accesible by any class
class RunParameters {
private:  
    
    RunParameters(); 
    ~RunParameters() {}; 

    /// The associated messenger
    UserMessenger<RunParameters> *fMessenger{nullptr}; 

    G4String f_outfile_path; 
    G4String f_infile_path; 
    
    ArmMode::EMode f_arm_mode;

    G4double f_vertical_raster_amplitude, f_horizontal_raster_amplitude; 

    G4double f_sieve_x_low, f_sieve_x_high, f_sieve_y_low, f_sieve_y_high;
    G4double f_momentum_low, f_momentum_high; 

    TargetMode::Bit f_target_mode{TargetMode::kNone}; 

    G4int f_expected_n_events_kept; 

public: 
    //_________________________________________________________________________
    // 
    //   List setters to set run information 
    //
    void Set_OutfilePath(G4String _x)   { f_outfile_path=_x; }
    void Set_InfilePath(G4String _x)    { f_infile_path=_x; }
    void Set_ArmMode(G4String); 
    void Set_SieveXLow(G4double _x)     { f_sieve_x_low=_x; }
    void Set_SieveXHigh(G4double _x)    { f_sieve_x_high=_x; }
    void Set_SieveYLow(G4double _x)     { f_sieve_y_low=_x; }
    void Set_SieveYHigh(G4double _x)    { f_sieve_x_high=_x; }    
    void Set_MomentumLow(G4double _x)   { f_momentum_low=_x; }
    void Set_MomentumHigh(G4double _x)  { f_momentum_high=_x; }
    void Set_VerticalRasterAmplitude(G4double _x)   { f_vertical_raster_amplitude=_x; }
    void Set_HorizontalRasterAmplitude(G4double _x) { f_horizontal_raster_amplitude=_x; }
    void Set_ExpectedNEventsKept(G4int _x)          { f_expected_n_events_kept=_x; }
    void Set_TargetMode(G4String); 


    G4String Get_OutfilePath() const    { return f_outfile_path; }  
    G4String Get_InfilePath() const     { return f_infile_path; }  
    ArmMode::EMode Get_ArmMode() const  { return f_arm_mode; } 
    G4double Get_SieveXLow() const      { return f_sieve_x_low; }
    G4double Get_SieveXHigh() const     { return f_sieve_x_high; }
    G4double Get_SieveYLow() const      { return f_sieve_y_low; }
    G4double Get_SieveYHigh() const     { return f_sieve_x_high; }    
    G4double Get_MomentumLow() const    { return f_momentum_low; }
    G4double Get_MomentumHigh() const   { return f_momentum_high; }
    G4double Get_VerticalRasterAmplitude() const    { return f_vertical_raster_amplitude; }
    G4double Get_HorizontalRasterAmplitude() const  { return f_horizontal_raster_amplitude; }
    G4int    Get_ExpectedNEventsKept() const        { return f_expected_n_events_kept; }
    TargetMode::Bit Get_TargetMode() const          { return f_target_mode; }
    
    bool RHRS_is_active() const { return f_arm_mode & ArmMode::kRHRS; }
    bool LHRS_is_active() const { return f_arm_mode & ArmMode::kLHRS; }
     
    //delete copy constructor
    RunParameters(const RunParameters&) = delete; 
    
    //delete copy assignment operator
    RunParameters& operator=(const RunParameters&) = delete; 


    /// @return ptr to singleton instance of global RunParameters object 
    static RunParameters& Instance()
    {
        static RunParameters inst; 
        return inst; 
    }; 

}; 

#endif 
