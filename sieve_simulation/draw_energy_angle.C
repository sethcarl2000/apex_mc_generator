#ifndef draw_energy_angle_C
#define draw_energy_angle_C

#include <TH2D.h>
#include <ROOT/RDataFrame.hxx>
#include <cmath> 
#include <string>
#include <functional> 
#include <stdexcept> 
#include <map> 

#include "TString.h"

/// @return a RDataFrame-ready filter that will filter particles based on charge
std::function<bool(double)> make_charge_filter(std::string particle_type);

/// @brief draw energy-angle distribution of tracks at the sieve
TH2D* draw_energy_angle(
    const char* path_infile,
    double cos_min = 0.988, 
    double cos_max = 1.000, 
    double E_min = 990.,
    double E_max = 2210., 
    std::string particle_type="all"
)
{
    /*auto h = new TH2D("h_energy_angle", "Energy - angle distribution of tracks at sieve;1 - cos(#theta);Energy (MeV)",
        200, 1 - cos_max, 1 - cos_min, 
        200, E_min, E_max
    );*/

    ROOT::RDataFrame df("tracks_sieve", path_infile); 

    auto h_cpy = df 
        
        .Filter(make_charge_filter(particle_type), {"charge"})

        .Define("energy", [](double px, double py, double pz){ return std::sqrt(px*px + py*py + pz*pz); },
        {"momentum_sieve_x", "momentum_sieve_y", "momentum_sieve_z"})
        
        .Define("one_minus_cos_theta", [](double px, double py, double pz)
        {
            return 1. - pz/std::sqrt(px*px + py*py + pz*pz);
        }, {"momentum_sieve_x", "momentum_sieve_y", "momentum_sieve_z"})

        .Histo2D<double>({"h_energy_angle_ptr", "Energy - angle distribution of tracks at sieve;1 - cos(#theta);Energy (MeV)",
            200, 1 - cos_max, 1 - cos_min, 
            200, E_min, E_max
        }, "one_minus_cos_theta", "energy"); 

    
    return (TH2D*)h_cpy->Clone("h_energy_angle");
}

std::function<bool(double)> make_charge_filter(std::string particle_type)
{
    if (particle_type=="electron") { return [](double charge){ return charge < 0.; }; }
    if (particle_type=="positron") { return [](double charge){ return charge > 0.; }; }
    if (particle_type=="photon")   { return [](double charge){ return charge==0.; }; } 

    if (particle_type=="all")      { return [](double charge){ return true; }; }

    throw std::logic_error("invalid particle type given");
    return [](double charge) { return false; };
}


#endif