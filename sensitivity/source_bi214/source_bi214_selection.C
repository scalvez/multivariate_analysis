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
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <stdlib.h>

void source_bi214_selection()
{

  // TFile * f_bi214 = TFile::Open("$SW_WORK_DIR/multivariate_analysis/data/.root");
  // TFile * f_bi214 = TFile::Open("bi214_tree.root");
  // TFile * f_bi214 = TFile::Open("pseudo.root");
  TFile * f_tl208 = TFile::Open("tl208_tree.root");

  // TTree *tree_bi214 = (TTree*)f_bi214->Get("snemodata");
  TTree *tree_tl208 = (TTree*)f_tl208->Get("snemodata");

  TCut cut_electron_energy = "1e_electron_energy > 1";
  // // TCut cut_electron_angle = "1e_electron_angle > 1.5";
  TCut cut_electron_angle = "2e_electrons_angle > 1.5";

  // tree_bi214->Draw("1e_electron_energy",cut_electron_energy + cut_electron_angle);
  TH1F* h = new TH1F("h","h",100,0,5);

  // tree_bi214->Project("h","1e_electron_energy",cut_electron_energy);
  // tree_bi214->Project("h","1e1a_alpha_track_length");
  // tree_bi214->Project("h","1e1a_alpha_track_length","","",1000);
  // tree_bi214->Project("h","1e1g_electron_gamma_energy_sum");
  // tree_bi214->Project("h","1e1g_electron_gamma_energy_sum","","",6760);
  // tree_tl208->Project("h","1e1g_electron_gamma_energy_sum");
  tree_tl208->Project("h","1e1g_electron_gamma_energy_sum","","",1535);

  TFile *f_output= new TFile("pseudo_tl208.root","RECREATE");

  // std::cout << h->Integral(0,100) << std::endl;

  // h->Scale(1./h->Integral(0,100));
  h->Write();
  return;
}
