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

#include "../config_sensitivity.h"
#include "spectra_1e.C"

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

void sensitivity_1e()
{
  // spectra_1e();

  TFile * f = TFile::Open("spectra_1e.root");

  TH1F *h_0nu = (TH1F*)f->Get("0nu");
  TH1F *h_2nu = (TH1F*)f->Get("2nu");
  TH1F *h_tl208 = (TH1F*)f->Get("tl208");
  TH1F *h_bi214 = (TH1F*)f->Get("bi214");
  TH1F *h_radon = (TH1F*)f->Get("radon");
  THStack *h_stack_norm = (THStack*)f->Get("hs");

  TGraphErrors * g_eff_0nu = new TGraphErrors();
  TGraphErrors * g_eff_2nu = new TGraphErrors();
  TGraphErrors * g_eff_tl208 = new TGraphErrors();
  TGraphErrors * g_eff_bi214 = new TGraphErrors();
  TGraphErrors * g_eff_radon = new TGraphErrors();

  TFile *f_output= new TFile("sensitivity_1e.root","RECREATE");

  TGraphErrors * significance = new TGraphErrors();
  significance->SetMarkerSize(2);
  significance->SetMarkerColor(kBlue);
  significance->SetDrawOption("AP");
  significance->SetTitle("Significance;Energy;Significance");

  TGraphErrors * halflife = new TGraphErrors();
  halflife->SetMarkerSize(2);
  halflife->SetMarkerColor(kBlue);
  halflife->SetDrawOption("AP");
  halflife->SetTitle("Halflife;Energy;Halflife");

  unsigned int nbins = h_0nu->GetNbinsX();

  double best_halflife_limit = 0;
  double Ncount = 0;

  for (unsigned int i = 0; i < nbins; ++i) {
    double S = 0;
    double B = 0;

    S =  h_0nu->Integral(i,nbins);
    B += h_2nu->Integral(i,nbins);
    B += h_tl208->Integral(i,nbins);
    B += h_bi214->Integral(i,nbins);
    B += h_radon->Integral(i,nbins);

    g_eff_0nu->SetPoint(i,h_0nu->GetBinLowEdge(i),h_0nu->Integral(i,nbins)*conf_sens::eff_0nu_1e);
    g_eff_2nu->SetPoint(i,h_2nu->GetBinLowEdge(i),h_2nu->Integral(i,nbins)*conf_sens::eff_2nu_1e);
    g_eff_tl208->SetPoint(i,h_tl208->GetBinLowEdge(i),h_tl208->Integral(i,nbins)*conf_sens::eff_tl208_1e);
    g_eff_bi214->SetPoint(i,h_bi214->GetBinLowEdge(i),h_bi214->Integral(i,nbins)*conf_sens::eff_bi214_1e);
    g_eff_radon->SetPoint(i,h_radon->GetBinLowEdge(i),h_radon->Integral(i,nbins)*conf_sens::eff_radon_1e);

    if(S+B != 0)
      significance->SetPoint(i,h_0nu->GetBinLowEdge(i),S/sqrt(S+B));
    else
      significance->SetPoint(i,h_0nu->GetBinLowEdge(i),0);

    double eff_0nu = h_0nu->Integral(i,nbins); // normalized to 1

    double eff_2nu = h_2nu->Integral(i,nbins); // normalized to 1
    double N_2nu = eff_2nu * conf_sens::eff_2nu_1e * conf_sens::k_sens / conf_sens::T_2nu;
    double eff_tl208 = h_tl208->Integral(i,nbins);
    double N_tl208 = eff_tl208 * conf_sens::eff_tl208_1e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::tl208_activity;
    double eff_bi214 = h_bi214->Integral(i,nbins);
    double N_bi214 = eff_bi214 * conf_sens::eff_bi214_1e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::bi214_activity;
    double eff_radon = h_radon->Integral(i,nbins);
    double N_radon = eff_radon * conf_sens::eff_radon_1e * conf_sens::tracker_volume * conf_sens::exposure * conf_sens::year2sec * conf_sens::radon_activity;

    double N_excluded = get_number_of_excluded_events(N_2nu + N_tl208 + N_bi214 + N_radon);
    // double N_excluded = get_number_of_excluded_events(N_2nu + N_tl208 + N_bi214);
    double tmp_halflife = eff_0nu * conf_sens::eff_0nu_1e * conf_sens::k_sens / N_excluded;
    if(best_halflife_limit<tmp_halflife) {
      Ncount = N_excluded;
    }

    best_halflife_limit = std::max(best_halflife_limit, tmp_halflife);
    halflife->SetPoint(i,h_0nu->GetBinLowEdge(i),tmp_halflife);
  }

  // significance->Draw("AP");
  significance->SetName("significance_norm");
  significance->Write();

  TDirectory *isotope = f_output->mkdir("isotope");
  isotope->cd();

  h_0nu->SetName("h_0nu_norm");
  h_2nu->SetName("h_2nu_norm");
  h_tl208->SetName("h_tl208_norm");
  h_bi214->SetName("h_bi214_norm");
  h_radon->SetName("h_radon_norm");
  h_0nu->Write();
  h_2nu->Write();
  h_tl208->Write();
  h_bi214->Write();
  h_radon->Write();

  THStack *hs_count = new THStack("hs_count","Stacked spectra");
  h_0nu->Scale(10);
  h_2nu->Scale(conf_sens::eff_2nu_1e * conf_sens::k_sens / conf_sens::T_2nu );
  h_tl208->Scale(conf_sens::eff_tl208_1e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::tl208_activity );
  h_bi214->Scale(conf_sens::eff_bi214_1e * conf_sens::isotope_mass * conf_sens::exposure * conf_sens::year2sec * conf_sens::bi214_activity );
  h_radon->Scale(conf_sens::eff_radon_1e * conf_sens::tracker_volume * conf_sens::exposure * conf_sens::year2sec * conf_sens::radon_activity );

  hs_count->Add(h_2nu);
  hs_count->Add(h_tl208);
  hs_count->Add(h_bi214);
  hs_count->Add(h_radon);
  hs_count->Add(h_0nu);

  f_output->cd();
  hs_count->Write();

  isotope->cd();

  h_0nu->SetName("h_0nu_count");
  h_2nu->SetName("h_2nu_count");
  h_tl208->SetName("h_tl208_count");
  h_bi214->SetName("h_bi214_count");
  h_radon->SetName("h_radon_count");
  h_0nu->Write();
  h_2nu->Write();
  h_tl208->Write();
  h_bi214->Write();
  h_radon->Write();

  std::cout << " Best halflife " << std::endl;
  std::cout << "     T_1/2 :  "  << best_halflife_limit << "  ,   N_bg = " << Ncount << std::endl;

  f_output->cd();
  halflife->SetName("halflife");
  halflife->Write();

  TDirectory *efficiencies = f_output->mkdir("efficiencies");
  efficiencies->cd();

  g_eff_0nu->SetName("g_eff_0nu");
  g_eff_0nu->SetTitle("0#nu efficiency;Energy;Efficiency");
  g_eff_0nu->SetDrawOption("AL");
  g_eff_0nu->SetLineColor(kRed);
  g_eff_0nu->SetLineWidth(2);
  g_eff_0nu->Write();

  g_eff_2nu->SetName("g_eff_2nu");
  g_eff_2nu->SetTitle("2#nu efficiency;Energy;Efficiency");
  g_eff_2nu->SetDrawOption("AL");
  g_eff_2nu->SetLineColor(kBlue);
  g_eff_2nu->SetLineWidth(2);
  g_eff_2nu->Write();

  g_eff_tl208->SetName("g_eff_tl208");
  g_eff_tl208->SetTitle("^{208}Tl efficiency;Energy;Efficiency");
  g_eff_tl208->SetDrawOption("AL");
  g_eff_tl208->SetLineColor(kGreen+1);
  g_eff_tl208->SetLineWidth(2);
  g_eff_tl208->Write();

  g_eff_bi214->SetName("g_eff_bi214");
  g_eff_bi214->SetTitle("^{214}Bi efficiency;Energy;Efficiency");
  g_eff_bi214->SetDrawOption("AL");
  g_eff_bi214->SetLineColor(kOrange-3);
  g_eff_bi214->SetLineWidth(2);
  g_eff_bi214->Write();

  g_eff_radon->SetName("g_eff_radon");
  g_eff_radon->SetTitle("Radon efficiency;Energy;Efficiency");
  g_eff_radon->SetDrawOption("AL");
  g_eff_radon->SetLineColor(kMagenta);
  g_eff_radon->SetLineWidth(2);
  g_eff_radon->Write();

  // gStyle->SetTitleFontSize(0.08);
  // gPad->SetLogy();
  // // gPad->SetPad(0.1,0.2,0.9,0.95);
  // gPad->SetBottomMargin(0.155);
  // gPad->SetTopMargin(1);
  // gPad->SetRightMargin(1000);

  return;
}
