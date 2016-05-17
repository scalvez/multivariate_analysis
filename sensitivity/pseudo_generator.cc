#include "TH1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom.h"
#include <vector>
#include <map>
#include <iostream>

#include "sensitivity_constants.h"
// #include "sensitivity_measurements.h"

//to do : issue with redefinition of quantity_efficiency
// #include "channel_selection.cc"

void pseudo_generator(TString isotope, std::vector<TString> quantities, double activity, std::map < TString , double > & quantity_efficiency) {
// void pseudo_generator(TString isotope, std::vector<TString> quantities, double activity) {

  TString input_file = "../" + isotope + "_pdf.root";
  TString output_file = "../" + isotope + "_pseudo.root";

  TFile *f = TFile::Open(input_file);
  TFile *f_output= new TFile(output_file, "RECREATE");

  for(auto it = quantities.begin(); it != quantities.end(); ++it) {
    TString qty = *it;
    TH1 *h = (TH1*)f->Get(qty);
    TH1 *h_cdf = h->GetCumulative();
    // h_cdf->Draw();
    // h->Draw("same");
    TString key = isotope + "_" + qty;

    int n_events = int(activity*exposure*mass*quantity_efficiency.at(key));
    std::cout << "efficiency " << quantity_efficiency.at(key) << std::endl;
    std::cout << "  events " << n_events << std::endl;
    // int n_events = 10000;
    int nbins = h_cdf->GetNbinsX();

    TH1F *h_pseudo = new TH1F(qty,qty,nbins,h->GetXaxis()->GetBinLowEdge(1),h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));

    std::cout << "Simulating pseudo-experiment with " << n_events << std::endl;

    for(int i=0; i<n_events;++i) {
      double rdm = gRandom->Uniform(1);
      for(int j=1; j<=nbins;++j) {
        if(h_cdf->GetBinContent(j)>rdm) {
          h_pseudo->Fill(h_cdf->GetXaxis()->GetBinCenter(j));
          break;
        }
      }
    }

    // h->Scale(n_events);
    // h->Draw();
    // h_pseudo->SetLineColor(kRed);
    // h_pseudo->Draw("sameEP");

    h_pseudo->Write();
  }
  f->Close();
  f_output->Close();
  return;
}
