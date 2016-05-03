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

void fcn_to_minimize(int& npar, double* deriv, double& f, double par[], int flag)
{
  TFile * f_bi214 = TFile::Open("pdf.root");

  TH1F *bi214_alpha_track_length_pdf = (TH1F*)f_bi214->Get("h");

  TFile * f_bi214_pseudo = TFile::Open("pseudo.root");

  TH1F *bi214_alpha_track_length_pseudo = (TH1F*)f_bi214->Get("h");

  double channel_1e1a_efficiency = 10733./2500000;

  f = 0.;
  for(int i = 0; i < 100; i++)
    {
      double b = par[0]*channel_1e1a_efficiency*bi214_alpha_track_length_pdf->GetBinContent(i)*7*3.14e7;
      double d = bi214_alpha_track_length_pseudo->GetBinContent(i);
      f += 2*(b-d*log(b+0.000000001)+log(factorial(d)));
    }

  // std::cout << " DEBUG " << f << std::endl;
  return;

}

void source_bi214_fit()
{
  const int npar = 1;

  TMinuit minuit(npar);

  minuit.SetFCN(fcn_to_minimize);  //set the fcn to minimize

  double par[npar];               // the start values
  double stepSize[npar];          // step sizes
  double minVal[npar];            // minimum bound on parameter
  double maxVal[npar];            // maximum bound on parameter
  string parName[npar];           // parameter name


  par[0]= 1.;
  stepSize[0] = 10.;
  minVal[0] = 0;
  maxVal[0] = 1000;
  parName[0] = "activity in uBq/kg";

  minuit.DefineParameter(0, parName[0].c_str(),
                         par[0], stepSize[0], minVal[0], maxVal[0]);

  minuit.Migrad();

  double activity = 0;
  double activity_err = 0;
  minuit.GetParameter(0,activity,activity_err);

  std::cout << " Measured activity is : " << activity << " +/- " << activity_err << " uBq/kg" << std::endl;

  return;
}
