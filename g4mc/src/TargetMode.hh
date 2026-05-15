#ifndef TargetMode_HH
#define TargetMode_HH

#include <string> 

/// @brief bitflag to indicate target operation mode 
namespace TargetMode {
    
    enum Bit : int {
        kNone       = 0, 
        kProduction = 1<<0,
        kV1         = 1<<1,
        kV2         = 1<<2, 
        kV3         = 1<<3
    }; 

    std::string GetName(TargetMode::Bit bit) {
        switch (bit) {
            case TargetMode::kNone : return "none";
            case TargetMode::kProduction : return "produciton";
            case TargetMode::kV1 : return "V1"; 
            case TargetMode::kV2 : return "V2"; 
            case TargetMode::kV3 : return "V3"; 
        }
    }
}

#endif