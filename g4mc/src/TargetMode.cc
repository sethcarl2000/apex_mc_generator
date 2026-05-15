#include "TargetMode.hh"

namespace TargetMode
{

std::string GetName(TargetMode::Bit bit) {
    switch (bit) {
        case TargetMode::kNone : return "none";
        case TargetMode::kProduction : return "produciton";
        case TargetMode::kV1 : return "V1"; 
        case TargetMode::kV2 : return "V2"; 
        case TargetMode::kV3 : return "V3";
        default : return "n/a"; 
    }
}

}