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
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <stdlib.h>

void roi()
{

  TFile * f_spectrum = TFile::Open("spectra.root");

  TH1F *h_0nu_roi = (TH1F*)f_spectrum->Get("energy/Se82.0nubb_source_strips_bulk_");
  TH1F *h_2nu_roi = (TH1F*)f_spectrum->Get("energy/Se82.2nubb-2MeV_source_strips_bulk_");
  TH1F *h_tl208_roi = (TH1F*)f_spectrum->Get("energy/Tl208_source_strips_bulk_");
  TH1F *h_bi214_roi = (TH1F*)f_spectrum->Get("energy/Bi214_Po214_source_strips_bulk_");
  TH1F *h_radon_roi = (TH1F*)f_spectrum->Get("energy/Bi214_Po214_field_wire_surface_");

  h_0nu_roi->Scale(1./h_0nu_roi->GetEntries());
  h_2nu_roi->Scale(1./h_2nu_roi->GetEntries());
  h_tl208_roi->Scale(1./h_tl208_roi->GetEntries());
  h_bi214_roi->Scale(1./h_bi214_roi->GetEntries());
  h_radon_roi->Scale(1./h_radon_roi->GetEntries());

  TFile *f_output= new TFile("roi_spectra.root","RECREATE");

  h_0nu_roi->SetLineColor(kRed);
  h_0nu_roi->SetLineWidth(2);
  h_0nu_roi->SetFillColor(kRed);
  h_0nu_roi->SetTitle("0#nu;ROI energy; Probability");
  h_0nu_roi->SetName("0nu");
  // h_0nu->Rebin();

  h_2nu_roi->SetLineColor(kBlue);
  h_2nu_roi->SetLineWidth(2);
  h_2nu_roi->SetFillColor(kBlue);
  h_2nu_roi->SetTitle("2#nu;ROI energy; Probability");
  h_2nu_roi->SetName("2nu");
  // h_2nu->Rebin();

  h_tl208_roi->SetLineColor(kGreen+1);
  h_tl208_roi->SetLineWidth(2);
  h_tl208_roi->SetFillColor(kGreen+1);
  h_tl208_roi->SetTitle("^{208}Tl;ROI energy; Probability");
  h_tl208_roi->SetName("tl208");
  // h_tl208->Rebin();

  h_bi214_roi->SetLineColor(kOrange-3);
  h_bi214_roi->SetLineWidth(2);
  h_bi214_roi->SetFillColor(kOrange-3);
  h_bi214_roi->SetTitle("^{214}Bi;ROI energy; Probability");
  h_bi214_roi->SetName("bi214");
  // h_bi214->Rebin();

  h_radon_roi->SetLineColor(kMagenta);
  h_radon_roi->SetLineWidth(2);
  h_radon_roi->SetFillColor(kMagenta);
  h_radon_roi->SetTitle("Radon;ROI energy; Probability");
  h_radon_roi->SetName("radon");
  // h_radon->Rebin();

  Double_t xl1=.7, yl1=0.775, xl2=0.9, yl2=0.9;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->AddEntry(h_0nu_roi,"0#nu");
  leg->AddEntry(h_2nu_roi,"2#nu");
  leg->AddEntry(h_tl208_roi,"^{208}Tl");
  leg->AddEntry(h_bi214_roi,"^{214}Bi");
  leg->AddEntry(h_radon_roi,"Radon");
  leg->SetFillColor(kWhite);

  // h_0nu_roi->Draw("");
  // h_2nu_roi->Draw("same");
  // h_tl208_roi->Draw("same");
  // h_bi214_roi->Draw("same");
  // h_radon_roi->Draw("same");
  // leg->Draw("same");

  h_0nu_roi->Write();
  h_2nu_roi->Write();
  h_tl208_roi->Write();
  h_bi214_roi->Write();
  h_radon_roi->Write();

  THStack *hs = new THStack("hs","Stacked spectra");
  hs->Add(h_2nu_roi);
  hs->Add(h_tl208_roi);
  hs->Add(h_bi214_roi);
  hs->Add(h_radon_roi);
  hs->Add(h_0nu_roi);

  hs->Draw();
  leg->Draw("same");

  hs->Write();

  // gStyle->SetTitleFontSize(0.08);
  // gPad->SetLogy();
  // // gPad->SetPad(0.1,0.2,0.9,0.95);
  // gPad->SetBottomMargin(0.155);
  // gPad->SetTopMargin(1);
  // gPad->SetRightMargin(1000);

  return;
}

/*
    // --- Example of simple scan

    TFile f("simple_scans.root");
    TGraphErrors *h_82 =(TGraphErrors*)f.Get("82");
    Double_t xl1=.05, yl1=0.75, xl2=xl1+.3, yl2=yl1+.125;
    TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
    leg->AddEntry(h0nu,"0#nu2#beta no cuts");
    leg->AddEntry(h2nu,"2#nu2#beta no cuts");


  */
