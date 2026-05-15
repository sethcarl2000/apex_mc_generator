#ifndef EventWriter_HH
#define EventWriter_HH

#include "TrackData_t.hh"
#include "ArmMode.hh"
#include "UserMessenger.hh"
#include <string> 
#include <vector>

//generalized class to write events to a TFile (implemented as a meyers singleton)
class EventWriter {
private: 
    //IMPORTANT: this assumes that events are 'paired' together properly; that, if events want to be saved in 
    // both-arm mode, that the Left-Arm TrackData_t must be added in the same order as the Right-Arm
    std::vector<TrackData_t> fData_singleArm; 

    //when data for both arms is kept, then keep this track data. 
    struct BothArmData_t { TrackData_t R, L; };

    std::vector<BothArmData_t> fData_both; 

    TrackData_t fTrackData_R, fTrackData_L; 

    EventWriter(); 
    ~EventWriter() {};

    //estimated max number of events. 
    // there will not necessarily be a crash if this estimate is exceeded, but it will affect performance. 
    long int fMax_events_estimate; 

    std::string fPath_outfile; 

    ArmMode::EMode fArm_mode{ArmMode::kNone}; 

    UserMessenger<EventWriter> fMessenger{nullptr}; 

    bool fIs_initialized{false};

public:

    //delete copy constructor & copy assignment operator
    EventWriter(const EventWriter&) = delete; 
    EventWriter& operator=(const EventWriter&) = delete; 

    //Get access to static instance of EventWriter 
    static EventWriter& Instance () {
        static EventWriter writer; 
        return writer; 
    }   

    /// @brief Initialize the data-output
    /// @param _mode arm mode to operate in 
    /// @param path_outfile path to output file 
    /// @param max_events_estimate estimate of max number of events to save. a crash will not neccessarily occur if this value is exceeded, but it will affect performance
    void Initialize(ArmMode::EMode _mode, std::string path_outfile, long int max_events_estimate); 


    /// @return 'false' if not yet initialized
    bool IsInitialized() const { return fIs_initialized; }

    /// @return the ArmMode::EMode bit
    ArmMode::EMode GetArmMode() const { return fArm_mode; }
 
    /// @return Ref to track data struct for right / left arm 
    TrackData_t& GetTrackData(bool is_RHRS) { return is_RHRS ? fTrackData_R : fTrackData_L; }

    /// @return 'true' if RHRS is active
    bool RHRS_active() const { return GetArmMode() & ArmMode::kRHRS; }
    /// @return 'true' if LHRS is active 
    bool LHRS_active() const { return GetArmMode() & ArmMode::kLHRS; }

    //write event for one arm 
    void SaveEvent();

    /// @brief Save data to an output ROOT file 
    void WriteOutputFile(); 
};


#endif