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

void bdt_score_0nu_2nu()
{
  TFile * f_0nu = TFile::Open("bdt_score_0nu.root");
  TFile * f_2nu = TFile::Open("bdt_score_2nu.root");

  TH1F *h_0nu_bdt = (TH1F*)f_0nu->Get("MVA_BDT");
  TH1F *h_2nu_bdt = (TH1F*)f_2nu->Get("MVA_BDT");

  TFile *f_output= new TFile("bdt_scores_0nu_2nu.root","RECREATE");

  h_0nu_bdt->Scale(1./h_0nu_bdt->GetEntries());
  h_0nu_bdt->SetLineColor(kRed);
  h_0nu_bdt->SetLineWidth(2);
  h_0nu_bdt->SetFillColor(kRed);
  h_0nu_bdt->SetTitle("0#nu;BDT score; Probability");
  h_0nu_bdt->SetName("0nu");
  // h_0nu->Rebin();

  h_2nu_bdt->Scale(1./h_2nu_bdt->GetEntries());
  h_2nu_bdt->SetLineColor(kBlue);
  h_2nu_bdt->SetLineWidth(2);
  h_2nu_bdt->SetFillColor(kBlue);
  h_2nu_bdt->SetTitle("2#nu;BDT score; Probability");
  h_2nu_bdt->SetName("2nu");
  // h_2nu->Rebin();

  Double_t xl1=.75, yl1=0.7, xl2=0.9, yl2=0.9;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->AddEntry(h_0nu_bdt,"0#nu");
  leg->AddEntry(h_2nu_bdt,"2#nu");
  leg->SetFillColor(kWhite);

  h_0nu_bdt->Write();
  h_2nu_bdt->Write();

  THStack *hs = new THStack("hs","BDT scores");
  hs->Add(h_2nu_bdt);
  hs->Add(h_0nu_bdt);

  hs->SetTitle("BDT;BDT score;Probability");
  hs->Write();

  TCanvas * c1 = new TCanvas();
  c1->cd();
  hs->Draw();
  leg->Draw("same");
  c1->Write();

  return;
}
