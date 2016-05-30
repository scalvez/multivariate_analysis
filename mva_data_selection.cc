#include "Riostream.h"
#include "TMath.h"
#include "TH1.h"
#include "THStack.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"
#include "TTree.h"
#include "TBranch.h"
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
#include <limits>

void mva_data_selection()
{
  TFile *oldfile = TFile::Open("/home/calvez/nemo/work_dir/multivariate_analysis/sensitivity/trees/radon_tree.root");
  TFile *newfile = new TFile("./data_radon.root","recreate");

  TTree *oldtree = (TTree*)oldfile->Get("snemodata");

  TTree *newtree = new TTree("snemodata","snemodata");

  oldtree->SetBranchStatus("*",0);
  oldtree->SetBranchStatus("2e_*",1);

  Double_t old_electron_minimal_energy;
  Double_t old_electron_maximal_energy;
  Double_t old_electrons_energy_difference;
  Double_t old_electrons_energy_sum;
  Double_t old_electrons_internal_probability;
  Double_t old_electrons_external_probability;
  Double_t old_electrons_vertices_probability;
  Double_t old_electrons_angle;
  Double_t old_electrons_cos_angle;
  Double_t old_electron_Emin_track_length;
  Double_t old_electron_Emax_track_length;

  oldtree->SetBranchAddress("2e_electron_minimal_energy",&old_electron_minimal_energy);
  oldtree->SetBranchAddress("2e_electron_maximal_energy",&old_electron_maximal_energy);
  oldtree->SetBranchAddress("2e_electrons_energy_difference",&old_electrons_energy_difference);
  oldtree->SetBranchAddress("2e_electrons_energy_sum",&old_electrons_energy_sum);
  oldtree->SetBranchAddress("2e_electrons_internal_probability",&old_electrons_internal_probability);
  oldtree->SetBranchAddress("2e_electrons_external_probability",&old_electrons_external_probability);
  oldtree->SetBranchAddress("2e_electrons_vertices_probability",&old_electrons_vertices_probability);
  oldtree->SetBranchAddress("2e_electrons_angle",&old_electrons_angle);
  oldtree->SetBranchAddress("2e_electrons_cos_angle",&old_electrons_cos_angle);
  oldtree->SetBranchAddress("2e_electron_Emin_track_length",&old_electron_Emin_track_length);
  oldtree->SetBranchAddress("2e_electron_Emax_track_length",&old_electron_Emax_track_length);

  Double_t electron_minimal_energy = 0;
  Double_t electron_maximal_energy = 0;
  Double_t electrons_energy_difference = 0;
  Double_t electrons_energy_sum = 0;
  Double_t electrons_internal_probability = 0;
  Double_t electrons_external_probability = 0;
  Double_t electrons_vertices_probability = 0;
  Double_t electrons_angle = 0;
  Double_t electrons_cos_angle = 0;
  Double_t electron_Emin_track_length = 0;
  Double_t electron_Emax_track_length = 0;

  newtree->Branch("2e_electron_minimal_energy",&electron_minimal_energy,"2e_electron_minimal_energy/D");
  newtree->Branch("2e_electron_maximal_energy",&electron_maximal_energy,"2e_electron_maximal_energy/D");
  newtree->Branch("2e_electrons_energy_difference",&electrons_energy_difference,"2e_electrons_energy_difference/D");
  newtree->Branch("2e_electrons_energy_sum",&electrons_energy_sum,"2e_electrons_energy_sum/D");
  newtree->Branch("2e_electrons_internal_probability",&electrons_internal_probability,"2e_electrons_internal_probability/D");
  newtree->Branch("2e_electrons_external_probability",&electrons_external_probability,"2e_electrons_external_probability/D");
  newtree->Branch("2e_electrons_vertices_probability",&electrons_vertices_probability,"2e_electrons_vertices_probability/D");
  newtree->Branch("2e_electrons_angle",&electrons_angle,"2e_electrons_angle/D");
  newtree->Branch("2e_electrons_cos_angle",&electrons_cos_angle,"2e_electrons_cos_angle/D");
  newtree->Branch("2e_electron_Emin_track_length",&electron_Emin_track_length,"2e_electron_Emin_track_length/D");
  newtree->Branch("2e_electron_Emax_track_length",&electron_Emax_track_length,"2e_electron_Emax_track_length/D");

  int nentries = oldtree->GetEntriesFast();
  int n_wanted = 1e6;
  int count = 0;
  for(int i = 0; i< nentries; ++i) {
    oldtree->GetEntry(i);

    if (i%100000 == 0)
      std::cout << " i " << i << std::endl;

    //don't keep overflow
    if(old_electron_minimal_energy != old_electron_minimal_energy)
      continue;

    electron_minimal_energy = old_electron_minimal_energy;
    electron_maximal_energy = old_electron_maximal_energy;
    electrons_energy_difference = old_electrons_energy_difference;
    electrons_energy_sum = old_electrons_energy_sum;
    electrons_internal_probability = old_electrons_internal_probability;
    electrons_external_probability = old_electrons_external_probability;
    electrons_vertices_probability = old_electrons_vertices_probability;
    electrons_angle = old_electrons_angle;
    electrons_cos_angle = old_electrons_cos_angle;
    electron_Emin_track_length = old_electron_Emin_track_length;
    electron_Emax_track_length = old_electron_Emax_track_length;

    TBranch * b_2e_electron_minimal_energy = newtree->GetBranch("2e_electron_minimal_energy");
    // std::cout << i << "  " << old_electron_minimal_energy << std::endl;
    // std::cout << "      " << electron_minimal_energy << std::endl;

    // if(electron_minimal_energy > 5)
    //   continue;

    // b_2e_electron_minimal_energy->Fill();
    newtree->Fill();
    count++;
    if(count == n_wanted)
      break;
  }

  newfile->Write();
  delete oldfile;
  delete newfile;


  return;
}
