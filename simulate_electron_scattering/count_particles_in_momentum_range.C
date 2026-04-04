#include <TGraph.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>  
#include <TH2D.h> 
#include <TAxis.h>
#include <TStyle.h>
#include <math.h>
#include <ROOT/RDataFrame.hxx>
#include <vector> 
#include <iostream> 

using namespace std; 

namespace {
    constexpr double mom_min = 1104*( 1 - 0.06 ); 
    constexpr double mom_max = 1104*( 1 + 0.06 ); 
    
    constexpr double stat_scale = 4e5; 

    constexpr double beam_energy = 2200; //MeV
}

double Count_particles_in_acceptance(ROOT::RDF::RNode df, double p_min, double p_max, int charge); 

int count_particles_in_momentum_range()
{
    const int min_thickness_um = 50; 
    const int max_thickness_um = 12500; 
    const int thickness_spacing = 50; 

    const int n_pts = (max_thickness_um - min_thickness_um)/thickness_spacing; 

    const int nbins_energy = 60; 
    const double min_energy = 1000; 
    const double max_energy = beam_energy; 

    auto hist_energy_e = new TH2D(
        "h_energye", "Energy of punch-through lepton (e-);Target Thickness (mm);Energy (MeV)", 
        n_pts, ((double)min_thickness_um-(thickness_spacing/2))/1e3, ((double)max_thickness_um+(thickness_spacing/2))/1e3, 
        nbins_energy, min_energy, max_energy
    );
    
    auto hist_energy_p = new TH2D(
        "h_energyp", "Energy of punch-through lepton (e+);Target Thickness (mm);Energy (MeV)", 
        n_pts, min_thickness_um-(thickness_spacing/2), max_thickness_um+(thickness_spacing/2), 
        nbins_energy, min_energy, max_energy
    );

    vector<double> pts_thickness; 
    vector<double> pts_electron; 
    vector<double> pts_positron; 

    for (int um=min_thickness_um; um<=max_thickness_um; um+=thickness_spacing) {

        const char* path_file = Form("data/test_target_%ium.root",um);
        ROOT::RDataFrame df("tracks_sieve", path_file); 

        cout << "Processing file '" << path_file << "'" << endl; 
        
        pts_thickness.emplace_back(((double)um)/1e3); 
        //pts_electron .emplace_back(Count_particles_in_acceptance(df, mom_min, mom_max, -1) / stat_scale ); 
        //pts_positron .emplace_back(Count_particles_in_acceptance(df, mom_min, mom_max, +1) / stat_scale );
    
        auto df_energy = df
            
            .Define("energy", [](double px, double py, double pz)
            { 
                return sqrt(px*px + py*py + pz*pz); 
            }, {"momentum_sieve_x", "momentum_sieve_y", "momentum_sieve_z"});
            
        auto stats_e = *df_energy.Filter([](int charge, double p)
            { return (charge==-1) && (p > mom_min) && (p < mom_max); }, {"charge", "energy"}).Count(); 

        auto stats_p = *df_energy.Filter([](int charge, double p)
            { return (charge==+1) && (p > mom_min) && (p < mom_max); }, {"charge", "energy"}).Count(); 

        pts_electron.emplace_back( ((double)stats_e) / stat_scale ); 
        pts_positron.emplace_back( ((double)stats_p) / stat_scale ); 

        auto henergy_e = (TH1D*)df_energy 

            .Filter([](int charge){ return charge==-1; }, {"charge"})

            .Histo1D<double>({"hee", "energy", nbins_energy, min_energy, max_energy}, "energy")->Clone("hee_temp"); 
        henergy_e->SetDirectory(0);
        
        auto henergy_p = (TH1D*)df_energy 

            .Filter([](int charge){ return charge==+1; }, {"charge"})

            .Histo1D<double>({"hpe", "energy", nbins_energy, min_energy, max_energy}, "energy")->Clone("hpe_temp"); 
        henergy_p->SetDirectory(0);

        auto xax = henergy_e->GetXaxis(); 
        for (int bx=1; bx<=xax->GetNbins(); bx++) {

           double energy = xax->GetBinCenter(bx); 
           
           hist_energy_e->Fill( ((double)um)/1.e3, energy, henergy_e->GetBinContent(bx) / stat_scale );
           hist_energy_p->Fill( ((double)um)/1.e3, energy, henergy_p->GetBinContent(bx) / stat_scale ); 
        }
        delete henergy_e;
        delete henergy_p;
    }
    
    new TCanvas; 
    gPad->SetLogz(1);
    gStyle->SetOptStat(0);
    hist_energy_e->Draw("colz");
    
    new TCanvas; 
    gPad->SetLogz(1);
    hist_energy_p->Draw("colz");
    
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


