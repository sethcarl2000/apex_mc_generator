#ifndef ArmMode_HH
#define ArmMode_HH

#include <string> 

namespace ArmMode {

    enum EMode : int {
        kNone = 0,
        kRHRS = 1 << 0,
        kLHRS = 1 << 1, 
        kBoth = (kRHRS | kLHRS)
    }; 

    constexpr bool kRHRS_bool = true;
    constexpr bool kLHRS_bool = false;
    
    /// @return std::string corresponding to the name of a particular arm mode, for error reporting purposes. 
    std::string Name(EMode mode) {
        switch (mode) {
            case EMode::kBoth : return "Both"; 
            case EMode::kRHRS : return "RHRS"; 
            case EMode::kLHRS : return "LHRS";
            default : return "none";
        }
    }
}

#endif