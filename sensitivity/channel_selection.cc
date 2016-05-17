#include <iostream>

#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <map>

#include "analysis_config.h"
#include "channel_selection.h"
// #include "sensitivity_measurements.h"

void channel_selection(TString isotope, std::vector<TString> quantities_pdf, std::map < TString , double > & quantity_efficiency, bool normalize)
{
  TString input_file = "../" + isotope + "_tree.root";
  TString output_file = "../" + isotope + "_pdf.root";

  TFile *f = TFile::Open(input_file);
  TTree *tree = (TTree*)f->Get("snemodata");
  TFile *f_output= new TFile(output_file, "RECREATE");

  double isotope_mc_size = get_isotope_mc_size(input_file);

  // TCut cut = get_channel_cut(channel);

  for(unsigned int j=0; j<quantities_pdf.size(); ++j) {

    int nbins;
    double xmin, xmax;
    nbins = 100;
    xmin = 0;
    xmax = 5;

    TString qty = quantities_pdf.at(j);

    get_histogram_options(qty, nbins, xmin, xmax);

    TH1F* h = new TH1F(qty,qty,nbins,xmin,xmax);

    tree->Project(qty,qty);
    // tree->Project("h","1e1g_electron_gamma_energy_sum","","",1000);

    TString isotope_quantity = isotope + "_" + qty;
    // std::cout << "inserting " << isotope_quantity  <<  "  " <<  h->Integral(1,h->GetXaxis()->GetNbins())/isotope_mc_size << std::endl;
    quantity_efficiency.insert(std::pair<TString,double>(isotope_quantity,h->Integral(1,h->GetXaxis()->GetNbins())/isotope_mc_size));

    if(normalize)
      h->Scale(1./h->Integral(1,h->GetXaxis()->GetNbins()));
    h->Write();
  }

  f->Close();
  f_output->Close();

    return;
  }
