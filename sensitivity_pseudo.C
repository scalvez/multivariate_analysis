#include "Riostream.h"
#include "TMath.h"
#include "TH1.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "config_sensitivity.h"

void sensitivity_pseudo(TString bdt_file, TString spectrum_file, TString halflife_data)
{
  TFile * f_bdt = TFile::Open(bdt_file);
  TFile * f_roi = TFile::Open(spectrum_file);

  TH1F *h_bdt = (TH1F*)f_bdt->Get("MVA_BDT");
  TH1F *h_roi = (TH1F*)f_roi->Get("energy_spectrum");

  //MM
  // double local_eff_0nu_bdt_window = 0.62631 * 0.25462;
  // double local_eff_0nu_roi_window = 0.637263 * 0.25462;

  //RHC
  // double local_eff_0nu_bdt_window = 0.63023 * 0.161204;
  // double local_eff_0nu_roi_window = 0.611544 * 0.161204;

  // //MM extra bkg
  // double local_eff_0nu_bdt_window = 0.61512 * 0.25462;
  // double local_eff_0nu_roi_window = 0.637263 * 0.25462;

  //RHC extra bkg
  double local_eff_0nu_bdt_window = 0.59461 * 0.161204;
  double local_eff_0nu_roi_window = 0.611544 * 0.161204;


  int nbins_bdt = h_bdt->GetNbinsX();
  double N_excluded_bdt = get_number_of_excluded_events(h_bdt->Integral(71,nbins_bdt));
  double halflife_bdt = local_eff_0nu_bdt_window * conf_sens::k_sens / N_excluded_bdt;

  double N_excluded_roi = get_number_of_excluded_events(h_roi->Integral(55,61));
  double halflife_roi = local_eff_0nu_roi_window * conf_sens::k_sens / N_excluded_bdt;

  std::ofstream ofs;
  ofs.open (halflife_data, std::ofstream::out | std::ofstream::app);
  ofs << halflife_bdt << " " << halflife_roi << std::endl;
  ofs.close();

  // std::cout << "BDT "  << halflife_bdt << std::endl;
  // std::cout << "ROI "  << halflife_roi << std::endl;

  return;
}
