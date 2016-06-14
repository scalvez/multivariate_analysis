#include "Riostream.h"
#include "TMath.h"
#include "TH1.h"
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

void pseudo_halflives_distrib()
{
  TH1F *h_hl_bdt = new TH1F("h_hl_bdt","h_hl_bdt",120,2.5,6.5);
  TH1F *h_hl_roi = new TH1F("h_hl_roi","h_hl_roi",120,2.5,6.5);
  // TH1F *h_hl_bdt = new TH1F("h_hl_bdt","h_hl_bdt",120,0.5,4.5);
  // TH1F *h_hl_roi = new TH1F("h_hl_roi","h_hl_roi",120,0.5,4.5);

  char channel[200];
  int  nevent;

  TFile *f= new TFile("halflives_distrib.root","RECREATE");
  TCanvas *c= new TCanvas();
  char file_name[200];

  sprintf(file_name,"./halflives_pseudo_mm.txt");

  h_hl_bdt->SetName("BDT_halflives_distribution");
  h_hl_roi->SetName("ROI_halflives_distribution");

  std::fstream data;
  data.open(file_name, std::fstream::in);

  double hl_bdt =0;
  double hl_roi =0;

  while(data >> hl_bdt >> hl_roi)
    {
      h_hl_bdt->Fill(hl_bdt/1e24);
      h_hl_roi->Fill(hl_roi/1e24);
    }

  data.close();

  h_hl_bdt->SetLineColor(kRed);
  h_hl_bdt->SetLineWidth(2);
  h_hl_roi->SetLineColor(kBlue);
  h_hl_roi->SetLineWidth(2);

  h_hl_roi->SetLineWidth(2);

  h_hl_bdt->SetTitle(";Halflife limit [10^{24} years];# of pseudo-experiments");

  h_hl_bdt->SetStats(0);

  h_hl_bdt->Draw("");
  h_hl_roi->Draw("same");

  h_hl_bdt->Write();
  h_hl_roi->Write();

  Double_t xl1=.75, yl1=0.7, xl2=0.9, yl2=0.9;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->AddEntry(h_hl_bdt,"BDT");
  leg->AddEntry(h_hl_roi,"E_{TOT}");

  leg->Draw("same");

  std::cout << "mean bdt " << h_hl_bdt->GetMean(1) << std::endl;
  std::cout << "mean roi " << h_hl_roi->GetMean(1) << std::endl;
  // // h_cd->GetXaxis()->SetTitle("Channels");
  // h_cd->GetXaxis()->SetTitleFont(62);
  // h_cd->GetXaxis()->SetTitleSize(0.05);
  // h_cd->GetXaxis()->SetTitleOffset(0.88);
  // h_cd->GetXaxis()->SetRangeUser(0.,30.);
  // h_cd->GetYaxis()->SetTitle("Efficiency");
  // h_cd->GetYaxis()->SetTitleFont(62);
  // h_cd->GetYaxis()->SetTitleSize(0.05);
  // h_cd->GetYaxis()->SetTitleOffset(0.9);
  // h_cd->GetYaxis()->SetTitleOffset(0.9);
  // h_cd->GetYaxis()->SetRangeUser(0.00001,1.);
  // h_cd->SetMarkerStyle(20);
  // h_cd->SetMarkerSize(1);
  // h_cd->SetMarkerColor(kBlue);
  // h_cd->SetLineColor(kBlue);
  // h_cd->SetFillColor(kBlue);
  // h_cd->SetDrawOption("A");

  // h_cd->Draw("PE");
  // h_cd->Write();

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
