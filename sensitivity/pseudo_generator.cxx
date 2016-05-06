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
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TRandom.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "sensitivity_constants.h"

// void pseudo_generator(TString pdf_file= "", TString distribution = "") {
void pseudo_generator() {

  TString pdf_file = "bi214_pdf.root";
  TString distribution = "h";
  // pdf_file="";
  TFile * f = TFile::Open(pdf_file);

  // distribution = "";

  TH1 *h = (TH1*)f->Get(distribution);
  // h->Draw();
  TH1 *h_cdf = h->GetCumulative();
  // h_cdf->Draw();
  // h->Draw("same");

  double activity = 100e-6;
  int n_events = int(activity*bi214_channel_1e1g_efficiency*exposure*mass);
  // int n_events = 10000;
  int nbins = h_cdf->GetNbinsX();

  TH1F *h_pseudo = new TH1F("pseudo","pseudo",nbins,h->GetXaxis()->GetBinLowEdge(1),h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));

  std::cout << "Simulating pseudo-experiment with " << n_events << std::endl;

  for(unsigned int i=0; i<n_events;++i) {
    double rdm = gRandom->Uniform(1);
    for(unsigned int j=1; j<=nbins;++j) {
      if(h_cdf->GetBinContent(j)>rdm) {
        h_pseudo->Fill(h_cdf->GetXaxis()->GetBinCenter(j));
        break;
      }
    }
  }

  h->Scale(n_events);
  h->Draw();
  h_pseudo->SetLineColor(kRed);
  h_pseudo->Draw("sameEP");

  TFile *f_output= new TFile("pseudo.root","RECREATE");
  h_pseudo->Write();

  return;
}
