#include "ArmMode.hh"

namespace ArmMode
{

std::string Name(EMode mode) {
    switch (mode) {
        case EMode::kBoth : return "Both"; 
        case EMode::kRHRS : return "RHRS"; 
        case EMode::kLHRS : return "LHRS";
        default : return "none";
    }
}

};