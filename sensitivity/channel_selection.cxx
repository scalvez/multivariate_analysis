#include "Riostream.h"
#include "TMath.h"
#include "TH1.h"
#include "THStack.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLine.h"
#include "TTree.h"
#include "TLegend.h"
#include "TCut.h"
#include "TFile.h"
#include "TStyle.h"
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "analysis_config.h"

void get_histogram_options(TString quantity, int & nbins, double & xmin, double & xmax) {

  // hardcoded for now, see if it is parametrized in the configuration file
  if(quantity.Contains("energy")) {
    nbins = 100;
    xmin = 0;
    xmax = 5;
  }
  else if(quantity.Contains("probability")){
    nbins = 100;
    xmin = 0;
    xmax = 1;
  }
  else if(quantity.Contains("angle")) {
    nbins = 100;
    xmin = -1;
    xmax = 1;
  }
  else if (quantity.Contains("track_length")) {
    nbins = 100;
    xmin = 0;
    xmax = 500;
  }
  else {
    //Also maybe the alpha delayed time
    //for now, dummy values
    nbins = 100;
    xmin = 0;
    xmax = 1000;
  }

  return;
}

// TODO: to move in the analysis config and get MC dataset size from conf
double get_isotope_mc_size(TString isotope) {
  if(isotope.Contains("tl208")) {
    return 1e8;
  }
  else if(isotope.Contains("bi214")){
    return 1e8;
  }
  else {
    //tmp dirty
    return 1;
  }
}

void channel_selection(std::vector <TString> input_files, std::vector<TString> output_files, std::vector<TString> quantities_pdf, std::map <TString,double> qty_eff, bool normalize = true)
{
  for(unsigned int i = 0; i < input_files.size(); ++i) {
    TFile *f = TFile::Open(input_files.at(i));

    TTree *tree = (TTree*)f->Get("snemodata");

      TFile *f_output= new TFile(output_files.at(i),"RECREATE");

      double isotope_mc_size = get_isotope_mc_size(input_files.at(i));
      // TCut cut = get_channel_cut(channel);
      // TCut cut_electron_energy = "1e_electron_energy > 1";
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
        // tree->Project("h","1e1g_electron_gamma_energy_sum");

        qty_eff.insert (std::pair<TString,double>(qty,h->Integral(1,h->GetXaxis()->GetNbins())/isotope_mc_size));

        if(normalize)
          h->Scale(1./h->Integral(1,h->GetXaxis()->GetNbins()));
        h->Write();
      }
      f->Close();
      f_output->Close();
    }
  return;
}
