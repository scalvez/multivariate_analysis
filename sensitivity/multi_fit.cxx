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
#include "TMinuit.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "stats_tools.h"
#include "sensitivity_constants.h"

void fcn_to_minimize(int& npar, double* deriv, double& f, double par[], int flag)
{
  TFile * f_bi214 = TFile::Open("bi214_pdf.root");
  TH1F *bi214_pdf = (TH1F*)f_bi214->Get("h;1");

  TFile * f_tl208 = TFile::Open("tl208_pdf.root");
  TH1F *tl208_pdf = (TH1F*)f_tl208->Get("h;1");

  TFile * f_pseudo = TFile::Open("pseudo.root");
  TH1F *pseudo = (TH1F*)f_pseudo->Get("pseudo");

  f = 0.;
  for(int i = 0; i < 100; i++)
    {
      double b = (par[0]*bi214_pdf->GetBinContent(i)*bi214_channel_1e1g_efficiency
                  +par[1]*tl208_pdf->GetBinContent(i)*tl208_channel_1e1g_efficiency)*7*3.14e7*2.5;
      double d = pseudo->GetBinContent(i);

      if(b==0)
        b=1e-6;
      // std::cout << "i b d " << i << "  " << b << "  " << d << std::endl;
      f += 2*(b-d*log(b)+log_factorial(d));
    }

  f_bi214->Close();
  f_tl208->Close();
  f_pseudo->Close();
  return;

}

// double likelihood(double activity)
// {
//   TFile * f_bi214 = TFile::Open("bi214_pdf.root");

//   TH1F *bi214_alpha_track_length_pdf = (TH1F*)f_bi214->Get("h");

//   TFile * f_bi214_pseudo = TFile::Open("pseudo.root");

//   TH1F *bi214_alpha_track_length_pseudo = (TH1F*)f_bi214_pseudo->Get("h");

//   double f = 0.;

//   for(int i = 0; i < 100; i++)
//     {
//       double b = activity*bi214_channel_1e1a_efficiency*bi214_alpha_track_length_pdf->GetBinContent(i)*7*3.14e7*2.5;
//       int d = bi214_alpha_track_length_pseudo->GetBinContent(i);
//       // std::cout << "bin " << i << "  :  b = " << b << "  d = " << d << std::endl;
//       // std::cout << "bin likelihood " << 2*(b-d*log(b+0.1)+log(factorial(d))) << std::endl;
//       f += 2*(b-d*log(b+0.1)+log_factorial(d));
//     }

//   // std::cout << " DEBUG " << activity*1e6 << "  " << f << std::endl;
//   return f;
// }

void multi_fit()
{
  const int npar = 2;

  TMinuit minuit(npar);

  minuit.SetFCN(fcn_to_minimize);  //set the fcn to minimize

  double par[npar];               // the start values
  double stepSize[npar];          // step sizes
  double minVal[npar];            // minimum bound on parameter
  double maxVal[npar];            // maximum bound on parameter
  string parName[npar];           // parameter name

  par[0]= 10./1e6;
  stepSize[0] = 1e-6;
  minVal[0] = 1e-9;
  maxVal[0] = 1e-3;
  // parName[0] = "activity in uBq/kg";
  parName[0] = "activity bi214";

  minuit.DefineParameter(0, parName[0].c_str(),
                         par[0], stepSize[0], minVal[0], maxVal[0]);

  par[1]= 2./1e6;
  stepSize[1] = 1e-6;
  minVal[1] = 1e-9;
  maxVal[1] = 1e-3;
  // parName[1] = "activity in uBq/kg";
  parName[1] = "activity tl208";

  minuit.DefineParameter(1, parName[1].c_str(),
                         par[1], stepSize[1], minVal[1], maxVal[1]);

  minuit.Migrad();

  double activity_bi214 = 0;
  double activity_bi214_err = 0;
  minuit.GetParameter(0,activity_bi214,activity_bi214_err);

  double activity_tl208 = 0;
  double activity_tl208_err = 0;
  minuit.GetParameter(1,activity_tl208,activity_tl208_err);

  std::cout << " Measured bi214 activity is : " << activity_bi214*1e6 << " +/- " << activity_bi214_err*1e6 << " uBq/kg" << std::endl;
  std::cout << " Measured tl208 activity is : " << activity_tl208*1e6 << " +/- " << activity_tl208_err*1e6 << " uBq/kg" << std::endl;

  // TGraph *g_likelihood = new TGraph(100);
  // double activity_start = 1e-5;
  // double activity_end = 1e-3;

  // for (unsigned int i = 0; i<100;++i) {
  //   double activity = activity_start + i*(activity_end-activity_start)/100.;
  //   double Q = likelihood(activity);
  //   g_likelihood->SetPoint(i,activity*1e6,Q);
  // }
  // g_likelihood->Draw();

  TFile * f_bi214 = TFile::Open("bi214_pdf.root");
  TH1F *bi214_pdf = (TH1F*)f_bi214->Get("h");

  bi214_pdf->SetLineColor(kOrange);
  bi214_pdf->SetFillColor(kOrange);

  TFile * f_tl208 = TFile::Open("tl208_pdf.root");
  TH1F *tl208_pdf = (TH1F*)f_tl208->Get("h");
  tl208_pdf->SetLineColor(kGreen+1);
  tl208_pdf->SetFillColor(kGreen+1);

  TFile * f_pseudo = TFile::Open("pseudo.root");
  TH1F *pseudo = (TH1F*)f_pseudo->Get("pseudo");
  pseudo->SetLineColor(kBlack);
  // pseudo->SetMarkerStyle(1);

  THStack *hs = new THStack("hs","Energy [keV]");
  bi214_pdf->Scale(activity_bi214*bi214_channel_1e1g_efficiency*mass*exposure);
  tl208_pdf->Scale(activity_tl208*tl208_channel_1e1g_efficiency*mass*exposure);

  hs->Add(bi214_pdf);
  hs->Add(tl208_pdf);

  TCanvas *c1 = new TCanvas("c1","example",600,700);

  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  // pad1->SetBottomMargin(0.05);
  pad1->SetTopMargin(0.03);
  pad1->Draw();
  pad1->cd();

  hs->DrawClone();

  pseudo->DrawClone("samePE");
  c1->cd();

  // TH1 *h_sum = new TH1F("h_sum","h_sum");

  bi214_pdf->Sumw2();
  bi214_pdf->Add(tl208_pdf);
  //Error computation to check
  // bi214_pdf->Sumw2();
  pseudo->Sumw2();

  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0.0);
  // pad2->SetBottomMargin(0.05);
  pad2->Draw();
  pad2->cd();

  bi214_pdf->SetStats(0);
  pseudo->SetStats(0);

  pseudo->Divide(bi214_pdf);

  pseudo->GetYaxis()->SetRangeUser(0.6,1.4);
  pseudo->Draw("ep");

  TLine *line = new TLine(0,1,5,1);
  // line->SetLineColor(kBlack);
  line->SetLineStyle(kDashed);
  line->Draw("same");

  c1->cd();

  return;
}
