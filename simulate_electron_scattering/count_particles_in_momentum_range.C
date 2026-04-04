#include <TGraph.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>  
#include <math.h>
#include <ROOT/RDataFrame.hxx>
#include <vector> 
#include <iostream> 

using namespace std; 

namespace {
    constexpr double mom_min = 1104*( 1 - 0.06 ); 
    constexpr double mom_max = 1104*( 1 + 0.06 ); 
    
    constexpr double stat_scale = 4e5; 
}

double Count_particles_in_acceptance(ROOT::RDF::RNode df, double p_min, double p_max, int charge); 

int count_particles_in_momentum_range()
{
    const int min_thickness_um = 50; 
    const int max_thickness_um = 12500; 
    const int thickness_spacing = 50; 

    vector<double> pts_thickness; 
    vector<double> pts_electron; 
    vector<double> pts_positron; 

    for (int um=min_thickness_um; um<=max_thickness_um; um+=thickness_spacing) {

        const char* path_file = Form("data/test_target_%ium.root",um);
        ROOT::RDataFrame df("tracks_sieve", path_file); 

        cout << "Processing file '" << path_file << "'" << endl; 
        
        pts_thickness.emplace_back(((double)um)/1e3); 
        pts_electron .emplace_back(Count_particles_in_acceptance(df, mom_min, mom_max, -1) / stat_scale ); 
        pts_positron .emplace_back(Count_particles_in_acceptance(df, mom_min, mom_max, +1) / stat_scale );
    }
    auto legend = new TLegend; 

    new TCanvas; 
    auto ge = new TGraph(pts_thickness.size(), pts_thickness.data(), pts_electron.data()); 
    ge->SetTitle("Number of leptons in acceptance per Beam e-;Target thickness (mm);Count / beam e-");
    ge->SetLineColor(kBlack); 
    ge->SetLineWidth(2);
    ge->SetMinimum(0.);
    legend->AddEntry(ge, "e-"); 
    ge->Draw(); 
    
    auto gp = new TGraph(pts_thickness.size(), pts_thickness.data(), pts_positron.data()); 
    gp->SetLineColor(kRed); 
    gp->SetLineWidth(2); 
    //gp->SetLineStyle(kDashed); 
    legend->AddEntry(gp, "e+"); 
    gp->Draw("SAME"); 

    legend->Draw(); 

    return 0; 
}

double Count_particles_in_acceptance(ROOT::RDF::RNode df, double p_min, double p_max, int charge) 
{
    return (double)*df 
        .Filter([charge](int chg){ return chg==charge; }, {"charge"})
        .Define("momentum", [](double px, double py, double pz){ return sqrt(px*px + py*py + pz*pz); }, {"momentum_sieve_x","momentum_sieve_y","momentum_sieve_z"})
        .Filter([p_min, p_max](double p){ return (p > p_min) && (p < p_max); }, {"momentum"})
        .Count(); 
}


