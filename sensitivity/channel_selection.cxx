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

void channel_selection(TString input_file, TString channel="", TString output_file="channel_selection_output.root", bool normalize = true)
{
  TFile *f = TFile::Open(input_file);

  TTree *tree = (TTree*)f->Get("snemodata");

  // TCut cut = get_channel_cut(channel);
  // TCut cut_electron_energy = "1e_electron_energy > 1";

  TH1F* h = new TH1F("h","h",100,0,5);

  // tree->Project("h","1e1g_electron_gamma_energy_sum","","",1000);
  tree->Project("h","1e1g_electron_gamma_energy_sum");

  TFile *f_output= new TFile(output_file,"RECREATE");

  if(normalize)
    h->Scale(1./h->Integral(1,h->GetXaxis()->GetNbins()));

  h->Write();
  return;
}
