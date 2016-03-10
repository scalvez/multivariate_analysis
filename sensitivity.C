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

#include "config_sensitivity.h"

double get_number_of_excluded_events(const double number_of_events_)
{
  double number_of_excluded_events = 0.0;
  if (number_of_events_ < 29.0)
    {
      double x = number_of_events_;
      number_of_excluded_events =
        2.5617 + 0.747661 * x - 0.0666176 * std::pow(x,2)
        + 0.00432457 * std::pow(x,3) - 0.000139343 * std::pow(x,4)
        + 1.71509e-06 * std::pow(x,5);
    }
  else
    {
      number_of_excluded_events = 1.64 * std::sqrt(number_of_events_);
    }
  return number_of_excluded_events;
}

void sensitivity()
{

  TFile * f_bdt = TFile::Open("bdt_scores.root");

  TH1F *h_0nu_bdt = (TH1F*)f_bdt->Get("0nu");
  TH1F *h_2nu_bdt = (TH1F*)f_bdt->Get("2nu");
  TH1F *h_tl208_bdt = (TH1F*)f_bdt->Get("tl208");
  TH1F *h_bi214_bdt = (TH1F*)f_bdt->Get("bi214");
  TH1F *h_radon_bdt = (TH1F*)f_bdt->Get("radon");
  THStack *h_stack_bdt_norm = (THStack*)f_bdt->Get("hs");

  // TFile * f_spectrum = TFile::Open("spectra.root");

  // TH1F *h_0nu_roi = (TH1F*)f_spectrum->Get("energy/Se82.0nubb_source_strips_bulk_");
  // TH1F *h_2nu_roi = (TH1F*)f_spectrum->Get("energy/Se82.2nubb-2MeV_source_strips_bulk_");
  // TH1F *h_tl208_roi = (TH1F*)f_spectrum->Get("energy/Tl208_source_strips_bulk_");
  // TH1F *h_bi214_roi = (TH1F*)f_spectrum->Get("energy/Bi214_Po214_source_strips_bulk_");
  // TH1F *h_radon_roi = (TH1F*)f_spectrum->Get("energy/Bi214_Po214_field_wire_surface_");

  TFile * f_roi = TFile::Open("roi_spectra.root");

  TH1F *h_0nu_roi = (TH1F*)f_roi->Get("0nu");
  TH1F *h_2nu_roi = (TH1F*)f_roi->Get("2nu");
  TH1F *h_tl208_roi = (TH1F*)f_roi->Get("tl208");
  TH1F *h_bi214_roi = (TH1F*)f_roi->Get("bi214");
  TH1F *h_radon_roi = (TH1F*)f_roi->Get("radon");
  THStack *h_stack_roi_norm = (THStack*)f_roi->Get("hs");

  TFile *f_output= new TFile("sensitivity.root","RECREATE");



  TGraphErrors * bdt_significance = new TGraphErrors();
  bdt_significance->SetMarkerSize(2);
  bdt_significance->SetMarkerColor(kBlue);
  bdt_significance->SetDrawOption("AP");
  bdt_significance->SetTitle("BDT;BDT score;Significance");

  TGraphErrors * bdt_halflife = new TGraphErrors();
  bdt_halflife->SetMarkerSize(2);
  bdt_halflife->SetMarkerColor(kBlue);
  bdt_halflife->SetDrawOption("AP");
  bdt_halflife->SetTitle("BDT;BDT score;Halflife");

  unsigned int nbins = h_0nu_bdt->GetNbinsX();

  double best_halflife_limit_bdt = 0;

  for (unsigned int i = 0; i < nbins; ++i) {
    double S = 0;
    double B = 0;

    S =  h_0nu_bdt->Integral(i,nbins);
    B += h_2nu_bdt->Integral(i,nbins);
    B += h_tl208_bdt->Integral(i,nbins);
    B += h_bi214_bdt->Integral(i,nbins);
    B += h_radon_bdt->Integral(i,nbins);

    if(S+B != 0)
      bdt_significance->SetPoint(i,h_0nu_bdt->GetBinLowEdge(i),S/sqrt(S+B));
    else
      bdt_significance->SetPoint(i,h_0nu_bdt->GetBinLowEdge(i),0);

    double eff_0nu = h_0nu_bdt->Integral(i,nbins); // normalized to 1

    double eff_2nu = h_2nu_bdt->Integral(i,nbins); // normalized to 1
    double N_2nu = eff_2nu * conf_sens::eff_2nu_2e * conf_sens::k_sens / conf_sens::T_2nu;
    double eff_tl208 = h_tl208_bdt->Integral(i,nbins);
    double N_tl208 = eff_tl208 * conf_sens::eff_tl208_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::tl208_activity;
    double eff_bi214 = h_bi214_bdt->Integral(i,nbins);
    double N_bi214 = eff_bi214 * conf_sens::eff_bi214_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::bi214_activity;
    double eff_radon = h_radon_bdt->Integral(i,nbins);
    double N_radon = eff_radon * conf_sens::eff_radon_2e * conf_sens::tracker_volume * conf_sens::exposure * conf_sens::year2sec * conf_sens::radon_activity;

    double N_excluded = get_number_of_excluded_events(N_2nu + N_tl208 + N_bi214 + N_radon);
    double halflife = eff_0nu * conf_sens::k_sens / N_excluded;

    best_halflife_limit_bdt = std::max(best_halflife_limit_bdt, halflife);

    bdt_halflife->SetPoint(i,h_0nu_bdt->GetBinLowEdge(i),halflife);
  }

  // h_0nu_roi->Scale(1./h_0nu_roi->GetEntries());
  // h_2nu_roi->Scale(1./h_2nu_roi->GetEntries());
  // h_tl208_roi->Scale(1./h_tl208_roi->GetEntries());
  // h_bi214_roi->Scale(1./h_bi214_roi->GetEntries());
  // h_radon_roi->Scale(1./h_radon_roi->GetEntries());

  TGraphErrors * roi_significance = new TGraphErrors();
  roi_significance->SetMarkerSize(2);
  roi_significance->SetMarkerColor(kRed);
  roi_significance->SetDrawOption("AP");
  roi_significance->SetTitle("ROI;Energy;Significance");

  TGraphErrors * roi_halflife = new TGraphErrors();
  roi_halflife->SetMarkerSize(2);
  roi_halflife->SetMarkerColor(kBlue);
  roi_halflife->SetDrawOption("AP");
  roi_halflife->SetTitle("ROI;ROI energy;Halflife");

  double best_halflife_limit_roi = 0;
  unsigned int nbins_roi = 0;

  for (unsigned int i = 0; i < h_0nu_roi->GetNbinsX(); ++i) {
    double S = 0;
    double B = 0;
    S  = h_0nu_roi->Integral(i,nbins_roi);
    B += h_2nu_roi->Integral(i,nbins_roi);
    B += h_tl208_roi->Integral(i,nbins_roi);
    B += h_bi214_roi->Integral(i,nbins_roi);
    B += h_radon_roi->Integral(i,nbins_roi);

    if(S+B != 0)
      roi_significance->SetPoint(i,h_0nu_roi->GetBinLowEdge(i),S/sqrt(S+B));
    else
      roi_significance->SetPoint(i,h_0nu_roi->GetBinLowEdge(i),0);

     double eff_0nu = h_0nu_roi->Integral(i,nbins); // normalized to 1

    double eff_2nu = h_2nu_roi->Integral(i,nbins); // normalized to 1
    double N_2nu = eff_2nu * conf_sens::eff_2nu_2e * conf_sens::k_sens / conf_sens::T_2nu;
    double eff_tl208 = h_tl208_roi->Integral(i,nbins);
    double N_tl208 = eff_tl208 * conf_sens::eff_tl208_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::tl208_activity;
    double eff_bi214 = h_bi214_roi->Integral(i,nbins);
    double N_bi214 = eff_bi214 * conf_sens::eff_bi214_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::bi214_activity;
    double eff_radon = h_radon_roi->Integral(i,nbins);
    double N_radon = eff_radon * conf_sens::eff_radon_2e * conf_sens::tracker_volume * conf_sens::exposure * conf_sens::year2sec * conf_sens::radon_activity;

    double N_excluded = get_number_of_excluded_events(N_2nu + N_tl208 + N_bi214 + N_radon);
    double halflife = eff_0nu * conf_sens::k_sens / N_excluded;

    best_halflife_limit_roi = std::max(best_halflife_limit_roi, halflife);

    roi_halflife->SetPoint(i,h_0nu_roi->GetBinLowEdge(i),halflife);
  }

  // bdt_significance->Draw("AP");
  bdt_significance->SetName("bdt_significance_norm");
  roi_significance->SetName("roi_significance_norm");
  bdt_significance->Write();
  roi_significance->Write();

  h_stack_bdt_norm->SetName("hs_bdt_norm");
  h_stack_roi_norm->SetName("hs_roi_norm");
  h_stack_bdt_norm->Write();
  h_stack_roi_norm->Write();

  THStack *hs_bdt_count = new THStack("hs_bdt_count","Stacked spectra");
  h_2nu_bdt->Scale(conf_sens::eff_2nu_2e * conf_sens::k_sens / conf_sens::T_2nu );
  h_tl208_bdt->Scale(conf_sens::eff_tl208_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::tl208_activity );
  h_bi214_bdt->Scale(conf_sens::eff_bi214_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::bi214_activity );
  h_radon_bdt->Scale(conf_sens::eff_radon_2e * conf_sens::tracker_volume * conf_sens::exposure * conf_sens::year2sec * conf_sens::radon_activity );

  hs_bdt_count->Add(h_2nu_bdt);
  hs_bdt_count->Add(h_tl208_bdt);
  hs_bdt_count->Add(h_bi214_bdt);
  hs_bdt_count->Add(h_radon_bdt);
  hs_bdt_count->Add(h_0nu_bdt);

  hs_bdt_count->Write();

  THStack *hs_roi_count = new THStack("hs_roi_count","Stacked spectra");
  h_2nu_roi->Scale(conf_sens::eff_2nu_2e * conf_sens::k_sens / conf_sens::T_2nu );
  h_tl208_roi->Scale(conf_sens::eff_tl208_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::tl208_activity );
  h_bi214_roi->Scale(conf_sens::eff_bi214_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::bi214_activity );
  h_radon_roi->Scale(conf_sens::eff_radon_2e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::radon_activity );

  hs_roi_count->Add(h_2nu_roi);
  hs_roi_count->Add(h_tl208_roi);
  hs_roi_count->Add(h_bi214_roi);
  hs_roi_count->Add(h_radon_roi);
  hs_roi_count->Add(h_0nu_roi);

  hs_roi_count->Write();

  std::cout << " Best halflives " << std::endl;
  std::cout << "     BDT :  "  << best_halflife_limit_bdt << std::endl;
  std::cout << "     ROI :  "  << best_halflife_limit_roi << std::endl;

  bdt_halflife->SetName("bdt_halflife");
  roi_halflife->SetName("roi_halflife");
  bdt_halflife->Write();
  roi_halflife->Write();
  // h_0nu_bdt->Draw("");
  // h_2nu_bdt->Draw("same");
  // h_tl208_bdt->Draw("same");
  // h_bi214_bdt->Draw("same");
  // h_radon_bdt->Draw("same");
  // leg->Draw("same");

  // gStyle->SetTitleFontSize(0.08);
  // gPad->SetLogy();
  // // gPad->SetPad(0.1,0.2,0.9,0.95);
  // gPad->SetBottomMargin(0.155);
  // gPad->SetTopMargin(1);
  // gPad->SetRightMargin(1000);

  return;
}
