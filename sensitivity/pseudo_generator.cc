#include "TH1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom.h"
#include <vector>
#include <map>
#include <iostream>

#include "sensitivity_constants.h"
// #include "sensitivity_measurements.h"

extern std::map < TString , double > quantity_efficiency;

void pseudo_generator(TString isotope, std::vector<TString> quantities, double activity) {

  std::cout << " ---- Isotope : " << isotope << std::endl;

  TString input_file = "../" + isotope + "_pdf.root";
  TString output_file = "../" + isotope + "_pseudo.root";

  TFile *f = TFile::Open(input_file);
  TFile *f_output= new TFile(output_file, "RECREATE");

  TRandom *rdm = new TRandom();
  int n_decays = int(activity*exposure*mass);
  int n_decays_rdm = rdm->Poisson(n_decays);

  for(auto it = quantities.begin(); it != quantities.end(); ++it) {
    TString qty = *it;
    TH1 *h = (TH1*)f->Get(qty);
    TH1 *h_cdf = h->GetCumulative();
    TString key = isotope + "_" + qty;

    std::cout << " ---------- Qty : " << qty << std::endl;

    int n_events = int(n_decays_rdm*quantity_efficiency.at(key));

    // std::cout << " Generating pseudo-experiment with " << n_decays << " decays" << std::endl;
    // std::cout << " Generating pseudo-experiment with " << n_decays_rdm << " rdm decays" << std::endl;
    std::cout << " Generating pseudo-experiment with " << n_events << " events" << std::endl;
    int nbins = h_cdf->GetNbinsX();

    TH1F *h_pseudo = new TH1F(qty,qty,nbins,h->GetXaxis()->GetBinLowEdge(1),h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));

    for(int i=0; i<n_events;++i) {
      double rdm = gRandom->Uniform(1);
      for(int j=1; j<=nbins;++j) {
        if(h_cdf->GetBinContent(j)>rdm) {
          h_pseudo->Fill(h_cdf->GetXaxis()->GetBinCenter(j));
          break;
        }
      }
    }

    h_pseudo->Write();
  }
  f->Close();
  f_output->Close();
  return;
}
