#include "Riostream.h"
#include "TMath.h"
#include "TTree.h"
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

void comparison_spectra_calibration()
{

  TFile * f_0nu_nc = TFile::Open("./new_calib.root");
  TFile * f_0nu_pc = TFile::Open("./poisson_calib.root");

  TTree *tree_0nu_nc = (TTree*)f_0nu_nc->Get("snemodata");
  TTree *tree_0nu_pc = (TTree*)f_0nu_pc->Get("snemodata");

  TH1F *h_nc_emin = new TH1F("nc_emin","nc_emin",80,0,4);
  TH1F *h_pc_emin = new TH1F("pc_emin","pc_emin",80,0,4);
  TH1F *h_nc_emax = new TH1F("nc_emax","nc_emax",80,0,4);
  TH1F *h_pc_emax = new TH1F("pc_emax","pc_emax",80,0,4);
  TH1F *h_nc_ediff = new TH1F("nc_ediff","nc_ediff",80,0,4);
  TH1F *h_pc_ediff = new TH1F("pc_ediff","pc_ediff",80,0,4);
  TH1F *h_nc_esum = new TH1F("nc_esum","nc_esum",80,0,4);
  TH1F *h_pc_esum = new TH1F("pc_esum","pc_esum",80,0,4);

  double nc_emin = 0;
  double pc_emin = 0;
  double nc_emax = 0;
  double pc_emax = 0;
  double nc_ediff = 0;
  double pc_ediff = 0;
  double nc_esum = 0;
  double pc_esum = 0;

  tree_0nu_nc->SetBranchAddress("2e_electron_minimal_energy",&nc_emin);
  tree_0nu_nc->SetBranchAddress("2e_electron_maximal_energy",&nc_emax);
  tree_0nu_nc->SetBranchAddress("2e_electrons_energy_difference",&nc_ediff);
  tree_0nu_nc->SetBranchAddress("2e_electrons_energy_sum",&nc_esum);
  tree_0nu_pc->SetBranchAddress("2e_electron_minimal_energy",&pc_emin);
  tree_0nu_pc->SetBranchAddress("2e_electron_maximal_energy",&pc_emax);
  tree_0nu_pc->SetBranchAddress("2e_electrons_energy_difference",&pc_ediff);
  tree_0nu_pc->SetBranchAddress("2e_electrons_energy_sum",&pc_esum);

  int nentries_nc = tree_0nu_nc->GetEntriesFast();
  for(int i = 0; i< 1000000; ++i) {
    tree_0nu_nc->GetEntry(i);
    h_nc_emin->Fill(nc_emin);
    h_nc_emax->Fill(nc_emax);
    h_nc_ediff->Fill(nc_ediff);
    h_nc_esum->Fill(nc_esum);
  }

  int nentries_pc = tree_0nu_pc->GetEntriesFast();
  for(int i = 0; i< 1000000; ++i) {
    tree_0nu_pc->GetEntry(i);
    h_pc_emin->Fill(pc_emin);
    h_pc_emax->Fill(pc_emax);
    h_pc_ediff->Fill(pc_ediff);
    h_pc_esum->Fill(pc_esum);
  }

  TFile *f_output= new TFile("comparison_spectra_calibration.root","RECREATE");

  // h_nc_emin->Scale(1./h_nc_emin->GetEntries());
  h_nc_emin->SetLineColor(kRed);
  h_nc_emin->SetLineWidth(2);
  // h_nc_emin->SetFillColor(kRed);
  h_nc_emin->SetTitle("Electron minimal energy;Energy [keV]; Probability");
  h_nc_emin->SetName("h_nc_emin");
  // h_0nu->Rebin();

  // h_nc_emax->Scale(1./h_nc_emax->GetEntries());
  h_nc_emax->SetLineColor(kRed);
  h_nc_emax->SetLineWidth(2);
  // h_nc_emax->SetFillColor(kRed);
  h_nc_emax->SetTitle("Electron maximal energy;Energy [keV]; Probability");
  h_nc_emax->SetName("h_nc_emax");
  // h_0nu->Rebin();

  // h_nc_ediff->Scale(1./h_nc_ediff->GetEntries());
  h_nc_ediff->SetLineColor(kRed);
  h_nc_ediff->SetLineWidth(2);
  // h_nc_ediff->SetFillColor(kRed);
  h_nc_ediff->SetTitle("Electrons energy difference;Energy [keV]; Probability");
  h_nc_ediff->SetName("h_nc_ediff");
  // h_0nu->Rebin();

  // h_nc_esum->Scale(1./h_nc_esum->GetEntries());
  h_nc_esum->SetLineColor(kRed);
  h_nc_esum->SetLineWidth(2);
  // h_nc_esum->SetFillColor(kRed);
  h_nc_esum->SetTitle("Electrons energy sum;Energy [keV]; Probability");
  h_nc_esum->SetName("h_nc_esum");
  // h_0nu->Rebin();


  // h_pc_emin->Scale(1./h_pc_emin->GetEntries());
  h_pc_emin->SetLineColor(kBlue);
  h_pc_emin->SetLineWidth(2);
  // h_pc_emin->SetFillColor(kBlue);
  h_pc_emin->SetTitle("Electron minimal energy;Energy [keV]; Probability");
  h_pc_emin->SetName("h_pc_emin");
  // h_0nu->Rebin();

  // h_pc_emax->Scale(1./h_pc_emax->GetEntries());
  h_pc_emax->SetLineColor(kBlue);
  h_pc_emax->SetLineWidth(2);
  // h_pc_emax->SetFillColor(kBlue);
  h_pc_emax->SetTitle("Electron maximal energy;Energy [keV]; Probability");
  h_pc_emax->SetName("h_pc_emax");
  // h_0nu->Rebin();

  // h_pc_ediff->Scale(1./h_pc_ediff->GetEntries());
  h_pc_ediff->SetLineColor(kBlue);
  h_pc_ediff->SetLineWidth(2);
  // h_pc_ediff->SetFillColor(kBlue);
  h_pc_ediff->SetTitle("Electrons energy difference;Energy [keV]; Probability");
  h_pc_ediff->SetName("h_pc_ediff");
  // h_0nu->Rebin();

  // h_pc_esum->Scale(1./h_pc_esum->GetEntries());
  h_pc_esum->SetLineColor(kBlue);
  h_pc_esum->SetLineWidth(2);
  // h_pc_esum->SetFillColor(kBlue);
  h_pc_esum->SetTitle("Electrons energy sum;Energy [keV]; Probability");
  h_pc_esum->SetName("h_pc_esum");
  // h_0nu->Rebin();

  Double_t xl1=.65, yl1=0.8, xl2=0.9, yl2=0.97;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->AddEntry(h_nc_emin,"Gaussian calibration");
  leg->AddEntry(h_pc_emin,"Poisson calibration");
  leg->SetFillColor(kWhite);

  TCanvas *c1 = new TCanvas("c1","example",600,700);

  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  // pad1->SetBottomMargin(0.05);
  pad1->SetTopMargin(0.03);
  pad1->Draw();
  pad1->cd();
  // h_pc_emin->DrawNormalized("",1);
  // h_nc_emin->DrawNormalized("same",1);
  // h_pc_emax->DrawNormalized("",1);
  // h_nc_emax->DrawNormalized("same",1);
  // h_pc_ediff->DrawNormalized("",1);
  // h_nc_ediff->DrawNormalized("same",1);
  h_pc_esum->DrawNormalized("",1);
  h_nc_esum->DrawNormalized("same",1);
  leg->Draw("same");
  c1->cd();

  h_nc_emin->Sumw2();
  h_pc_emin->Sumw2();
  h_nc_emax->Sumw2();
  h_pc_emax->Sumw2();
  h_nc_ediff->Sumw2();
  h_pc_ediff->Sumw2();
  h_nc_esum->Sumw2();
  h_pc_esum->Sumw2();

  h_nc_emin->Scale(1./h_nc_emin->GetEntries());
  h_nc_emax->Scale(1./h_nc_emax->GetEntries());
  h_nc_ediff->Scale(1./h_nc_ediff->GetEntries());
  h_nc_esum->Scale(1./h_nc_esum->GetEntries());

  h_pc_emin->Scale(1./h_pc_emin->GetEntries());
  h_pc_emax->Scale(1./h_pc_emax->GetEntries());
  h_pc_ediff->Scale(1./h_pc_ediff->GetEntries());
  h_pc_esum->Scale(1./h_pc_esum->GetEntries());


  h_nc_emin->Write();
  h_nc_emax->Write();
  h_nc_ediff->Write();
  h_nc_esum->Write();

  h_pc_emin->Write();
  h_pc_emax->Write();
  h_pc_ediff->Write();
  h_pc_esum->Write();

  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0.0);
  // pad2->SetBottomMargin(0.05);
  pad2->Draw();
  pad2->cd();

  // h_nc_emin->SetStats(0);
  // h_nc_emin->Sumw2();

  h_pc_emin->SetStats(0);
  h_pc_emin->Sumw2();
  h_pc_emax->SetStats(0);
  h_pc_emax->Sumw2();
  h_pc_ediff->SetStats(0);
  h_pc_ediff->Sumw2();
  h_pc_esum->SetStats(0);
  h_pc_esum->Sumw2();

  h_pc_emin->Divide(h_nc_emin);
  // h_pc_emin->SetLineColor(kBlack);
  h_pc_emax->Divide(h_nc_emax);
  // h_pc_emax->SetLineColor(kBlack);
  h_pc_ediff->Divide(h_nc_ediff);
  h_pc_esum->Divide(h_nc_esum);

  h_pc_emin->GetYaxis()->SetRangeUser(0.9,1.1);
  h_pc_emax->GetYaxis()->SetRangeUser(0.9,1.1);
  h_pc_ediff->GetYaxis()->SetRangeUser(0.9,1.1);
  h_pc_esum->GetYaxis()->SetRangeUser(0.9,1.1);

  h_pc_emin->GetYaxis()->SetTitle("#frac{#color[4]{Poisson}}{#color[2]{Gauss}}  ");
  h_pc_emin->GetYaxis()->SetTitleSize(0.07);
  h_pc_emin->GetYaxis()->SetTitleOffset(0.45);
  h_pc_emax->GetYaxis()->SetTitle("#frac{#color[4]{Poisson}}{#color[2]{Gauss}}  ");
  h_pc_emax->GetYaxis()->SetTitleSize(0.07);
  h_pc_emax->GetYaxis()->SetTitleOffset(0.45);
  h_pc_ediff->GetYaxis()->SetTitle("#frac{#color[4]{Poisson}}{#color[2]{Gauss}}  ");
  h_pc_ediff->GetYaxis()->SetTitleSize(0.07);
  h_pc_ediff->GetYaxis()->SetTitleOffset(0.45);
  h_pc_esum->GetYaxis()->SetTitle("#frac{#color[4]{Poisson}}{#color[2]{Gauss}}  ");
  h_pc_esum->GetYaxis()->SetTitleSize(0.07);
  h_pc_esum->GetYaxis()->SetTitleOffset(0.45);

  // h_pc_emin->Draw("ep");
  // h_pc_emax->Draw("ep");
  // h_pc_ediff->Draw("ep");
  h_pc_esum->Draw("ep");

  // h_pc_emin->SetLineColor(kBlue);

  TLine *line = new TLine(0,1,4,1);
  // line->SetLineColor(kBlack);
  line->SetLineStyle(kDashed);
  line->Draw("same");

  c1->cd();

  return;
}
