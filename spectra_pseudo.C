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

void spectra_pseudo(TString input, TString output)
{
  TFile *f = TFile::Open(input);

  TTree *tree = (TTree*)f->Get("snemodata");

  TH1F *h = new TH1F("energy_spectrum","enegry_spectrum",100,0,5);

  TFile *f_output= new TFile(output,"RECREATE");

  double electrons_energy_sum = 0;
  tree->SetBranchAddress("2e_electrons_energy_sum",&electrons_energy_sum);

  int nentries = tree->GetEntriesFast();
  for(int i = 0; i< nentries; ++i) {
    tree->GetEntry(i);
    h->Fill(electrons_energy_sum);
  }

  h->Write();

  return;
}
