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

void activity_accuracy_exposure () {

  TGraph *g = new TGraph(4);

  //2nu
  // g->SetPoint(0,1./365.25,1e-1);
  // g->SetPoint(0,1./12,1e-1);
  // g->SetPoint(0,1.,5.9e-2);
  // // g->SetPoint(0,2.5,);

  // //tl
  // g->SetPoint(0,1./12,40);
  // g->SetPoint(1,6./12,11.9);
  // g->SetPoint(2,1,8);
  // g->SetPoint(3,2.5,5.5);

  // //bi
  g->SetPoint(0,1./12,35.8);
  g->SetPoint(1,6./12,5.8);
  g->SetPoint(2,1,3.8);
  g->SetPoint(3,2.5,2.4);


  g->GetXaxis()->SetTitleFont(62);
  g->GetXaxis()->SetTitleSize(0.05);
  g->GetXaxis()->SetTitleOffset(0.88);
  g->GetXaxis()->SetTitle("Time (years)");
  // g->GetYaxis()->SetTitle("Relative uncertainty on T_{1/2}^{2#nu} in %");
  // g->GetYaxis()->SetTitle("Relative uncertainty on A(^{208}Tl) in %");
  g->GetYaxis()->SetTitle("Relative uncertainty on A(^{214}Bi) in %");
  g->GetYaxis()->SetTitleFont(62);
  g->GetYaxis()->SetTitleSize(0.05);
  g->GetYaxis()->SetTitleOffset(0.9);
  g->GetYaxis()->SetTitleOffset(0.9);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1);
  // g->SetMarkerColor(kBlue);
  // g->SetLineColor(kBlue);
  // g->SetFillColor(kBlue);
  g->SetMarkerColor(kOrange-3);
  g->SetLineColor(kOrange-3);
  g->SetFillColor(kOrange-3);
  g->SetDrawOption("A");

  g->Draw("");
  return;
}
