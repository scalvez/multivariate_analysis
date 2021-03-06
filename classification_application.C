#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;

void classification_application( TString myMethodList = "" )
{
#ifdef __CINT__
  gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

  //---------------------------------------------------------------

  // This loads the library
  TMVA::Tools::Instance();

  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;

  // --- Cut optimisation
  Use["Cuts"]            = 0;
  Use["CutsD"]           = 0;
  Use["CutsPCA"]         = 0;
  Use["CutsGA"]          = 0;
  Use["CutsSA"]          = 0;
  //
  // --- 1-dimensional likelihood ("naive Bayes estimator")
  Use["Likelihood"]      = 0;
  Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
  Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
  Use["LikelihoodKDE"]   = 0;
  Use["LikelihoodMIX"]   = 0;
  //
  // --- Mutidimensional likelihood and Nearest-Neighbour methods
  Use["PDERS"]           = 0;
  Use["PDERSD"]          = 0;
  Use["PDERSPCA"]        = 0;
  Use["PDEFoam"]         = 0;
  Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
  Use["KNN"]             = 0; // k-nearest neighbour method
  //
  // --- Linear Discriminant Analysis
  Use["LD"]              = 0; // Linear Discriminant identical to Fisher
  Use["Fisher"]          = 0;
  Use["FisherG"]         = 0;
  Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
  Use["HMatrix"]         = 0;
  //
  // --- Function Discriminant analysis
  Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
  Use["FDA_SA"]          = 0;
  Use["FDA_MC"]          = 0;
  Use["FDA_MT"]          = 0;
  Use["FDA_GAMT"]        = 0;
  Use["FDA_MCMT"]        = 0;
  //
  // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
  Use["MLP"]             = 0; // Recommended ANN
  Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
  Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
  Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
  Use["TMlpANN"]         = 0; // ROOT's own ANN
  //
  // --- Support Vector Machine
  Use["SVM"]             = 0;
  //
  // --- Boosted Decision Trees
  Use["BDT"]             = 1; // uses Adaptive Boost
  Use["BDTG"]            = 0; // uses Gradient Boost
  Use["BDTB"]            = 0; // uses Bagging
  Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
  //
  // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
  Use["RuleFit"]         = 0;
  // ---------------------------------------------------------------
  Use["Plugin"]          = 0;
  Use["Category"]        = 0;
  Use["SVM_Gauss"]       = 0;
  Use["SVM_Poly"]        = 0;
  Use["SVM_Lin"]         = 0;

  std::cout << std::endl;
  std::cout << "==> Start classification application" << std::endl;

  // Select methods (don't look at this code - not of interest)
  if (myMethodList != "") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

    std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
    for (UInt_t i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);

      if (Use.find(regMethod) == Use.end()) {
        std::cout << "Method \"" << regMethod
                  << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
          std::cout << it->first << " ";
        }
        std::cout << std::endl;
        return;
      }
      Use[regMethod] = 1;
    }
  }

  // --------------------------------------------------------------------------------------------------

  // --- Create the Reader object

  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

  // Create a set of variables and declare them to the reader
  // - the variable names MUST corresponds in name and type to those given in the weight file(s) used

  Float_t electron_minimal_energy = 0;
  Float_t electron_maximal_energy = 0;
  Float_t electrons_energy_sum = 0;
  Float_t electrons_energy_difference = 0;
  Float_t electrons_internal_probability = 0;
  Float_t electrons_external_probability = 0;
  Float_t electrons_vertices_probability = 0;
  Float_t electrons_cos_angle = 0;
  Float_t electron_Emin_track_length = 0;
  Float_t electron_Emax_track_length = 0;

  reader->AddVariable( "2e_electron_minimal_energy", & electron_minimal_energy );
  reader->AddVariable( "2e_electron_maximal_energy", & electron_maximal_energy );
  reader->AddVariable( "2e_electrons_energy_sum", &electrons_energy_sum );
  reader->AddVariable( "2e_electrons_energy_difference", &electrons_energy_difference );
  reader->AddVariable( "2e_electrons_internal_probability", &electrons_internal_probability );
  reader->AddVariable( "2e_electrons_external_probability", &electrons_external_probability );
  reader->AddVariable( "2e_electrons_vertices_probability", &electrons_vertices_probability );
  reader->AddVariable( "2e_electrons_cos_angle", &electrons_cos_angle );
  reader->AddVariable( "2e_electron_Emin_track_length", &electron_Emin_track_length );
  reader->AddVariable( "2e_electron_Emax_track_length", &electron_Emax_track_length );

  // --- Book the MVA methods

  TString dir    = "weights/";
  TString prefix = "classification";

  // Book method(s)
  for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
    if (it->second) {
      TString methodName = TString(it->first) + TString(" method");
      TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
      reader->BookMVA( methodName, weightfile );
    }
  }

  // Book output histograms
  UInt_t nbin = 100;
  TH1F   *histLk(0), *histLkD(0), *histLkPCA(0), *histLkKDE(0), *histLkMIX(0), *histPD(0), *histPDD(0);
  TH1F   *histPDPCA(0), *histPDEFoam(0), *histPDEFoamErr(0), *histPDEFoamSig(0), *histKNN(0), *histHm(0);
  TH1F   *histFi(0), *histFiG(0), *histFiB(0), *histLD(0), *histNn(0),*histNnbfgs(0),*histNnbnn(0);
  TH1F   *histNnC(0), *histNnT(0), *histBdt(0), *histBdtG(0), *histBdtD(0), *histRf(0), *histSVMG(0);
  TH1F   *histSVMP(0), *histSVML(0), *histFDAMT(0), *histFDAGA(0), *histCat(0), *histPBdt(0);

  if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
  if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
  if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
  if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
  if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
  if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
  if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
  if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
  if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
  if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
  if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
  if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
  if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
  if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
  if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
  if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
  if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
  if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
  if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
  if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -1.1, 1.1 );
  // if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
  if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
  if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
  if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
  if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
  if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
  if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
  if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
  if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
  if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
  if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

  // PDEFoam also returns per-event error, fill in histogram, and also fill significance
  if (Use["PDEFoam"]) {
    histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
    histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
    histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
  }

  // Book example histogram for probability (the other methods are done similarly)
  TH1F *probHistFi(0), *rarityHistFi(0);
  if (Use["Fisher"]) {
    probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
    rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
  }

  // Prepare input tree (this must be replaced by your data source)
  // in this example, there is a toy tree with signal and one with background events
  // we'll later on use only the "signal" events for the test in this example.
  //
TFile *input(0);

TString fname = "./data_radon.root";

input = TFile::Open( fname ); // check if file in local directory exists

  if (!input) {
    std::cout << "ERROR: could not open data file" << std::endl;
    exit(1);
  }
  std::cout << "--- Classification Application   : Using input file: " << input->GetName() << std::endl;

  // --- Event loop
  std::cout << "--- Select signal sample" << std::endl;
  TTree* theTree = (TTree*)input->Get("snemodata");

  Double_t br_2e_electron_minimal_energy = 1;
  Double_t br_2e_electron_maximal_energy = 1;
  Double_t br_2e_electrons_energy_difference = 1;
  Double_t br_2e_electrons_energy_sum = 1;
  Double_t br_2e_electrons_internal_probability = 1;
  Double_t br_2e_electrons_external_probability = 1;
  Double_t br_2e_electrons_vertices_probability = 1;
  Double_t br_2e_electrons_cos_angle = 1;
  Double_t br_2e_electron_Emin_track_length = 1;
  Double_t br_2e_electron_Emax_track_length = 1;

  theTree->SetBranchAddress( "2e_electron_minimal_energy", &br_2e_electron_minimal_energy);
  theTree->SetBranchAddress( "2e_electron_maximal_energy", &br_2e_electron_maximal_energy);
  theTree->SetBranchAddress( "2e_electrons_energy_difference", &br_2e_electrons_energy_difference);
  theTree->SetBranchAddress( "2e_electrons_energy_sum", &br_2e_electrons_energy_sum);
  theTree->SetBranchAddress( "2e_electrons_internal_probability", &br_2e_electrons_internal_probability);
  theTree->SetBranchAddress( "2e_electrons_external_probability", &br_2e_electrons_external_probability);
  theTree->SetBranchAddress( "2e_electrons_vertices_probability", &br_2e_electrons_vertices_probability);
  theTree->SetBranchAddress( "2e_electrons_cos_angle", &br_2e_electrons_cos_angle);
  theTree->SetBranchAddress( "2e_electron_Emin_track_length", &br_2e_electron_Emin_track_length);
  theTree->SetBranchAddress( "2e_electron_Emax_track_length", &br_2e_electron_Emax_track_length);

  // Efficiency calculator for cut method
  Int_t    nSelCutsGA = 0;
  Double_t effS       = 0.7;

  std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

  std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
  TStopwatch sw;
  sw.Start();

for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

    if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

    theTree->GetEntry(ievt);

    electron_minimal_energy = (Float_t) br_2e_electron_minimal_energy;
    electron_maximal_energy = (Float_t) br_2e_electron_minimal_energy;
    electrons_energy_difference = (Float_t) br_2e_electrons_energy_difference;
    electrons_energy_sum = (Float_t) br_2e_electrons_energy_sum;
    electrons_internal_probability = (Float_t) br_2e_electrons_internal_probability;
    electrons_external_probability = (Float_t) br_2e_electrons_external_probability;
    electrons_vertices_probability = (Float_t) br_2e_electrons_vertices_probability;
    electrons_cos_angle = (Float_t) br_2e_electrons_cos_angle;
    electron_Emin_track_length = (Float_t) br_2e_electron_Emin_track_length;
    electron_Emax_track_length = (Float_t) br_2e_electron_Emax_track_length;

    if (Use["CutsGA"]) {
      // Cuts is a special case: give the desired signal efficienciy
      Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
      if (passed) nSelCutsGA++;
    }

    if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
    if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
    if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
    if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
    if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
    if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
    if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
    if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
    if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
    if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
    if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
    if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
    if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
    if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
    if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
    if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
    if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
    if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
    if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
    if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
    if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
    if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) );
    if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
    if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
    if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
    if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
    if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
    if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
    if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
    if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );

    // Retrieve also per-event error
    if (Use["PDEFoam"]) {
      Double_t val = reader->EvaluateMVA( "PDEFoam method" );
      Double_t err = reader->GetMVAError();
      histPDEFoam   ->Fill( val );
      histPDEFoamErr->Fill( err );
      if (err>1.e-50) histPDEFoamSig->Fill( val/err );
    }

    // Retrieve probability instead of MVA output
    if (Use["Fisher"])   {
      probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) );
      rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) );
    }
  }

  // Get elapsed time
  sw.Stop();
  std::cout << "--- End of event loop: "; sw.Print();

  // Get efficiency for cuts classifier
  if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                               << " (for a required signal efficiency of " << effS << ")" << std::endl;

  if (Use["CutsGA"]) {

    // test: retrieve cuts for particular signal efficiency
    // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer
    TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

    if (mcuts) {
      std::vector<Double_t> cutsMin;
      std::vector<Double_t> cutsMax;
      mcuts->GetCuts( 0.7, cutsMin, cutsMax );
      std::cout << "--- -------------------------------------------------------------" << std::endl;
      std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
      for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
        std::cout << "... Cut: "
                  << cutsMin[ivar]
                  << " < \""
                  << mcuts->GetInputVar(ivar)
                  << "\" <= "
                  << cutsMax[ivar] << std::endl;
      }
      std::cout << "--- -------------------------------------------------------------" << std::endl;
    }
  }

  // --- Write histograms

TFile *target = new TFile ("bdt_score_radon_merge_3.root","RECREATE");

if (Use["Likelihood"   ])   histLk     ->Write();
  if (Use["LikelihoodD"  ])   histLkD    ->Write();
  if (Use["LikelihoodPCA"])   histLkPCA  ->Write();
  if (Use["LikelihoodKDE"])   histLkKDE  ->Write();
  if (Use["LikelihoodMIX"])   histLkMIX  ->Write();
  if (Use["PDERS"        ])   histPD     ->Write();
  if (Use["PDERSD"       ])   histPDD    ->Write();
  if (Use["PDERSPCA"     ])   histPDPCA  ->Write();
  if (Use["KNN"          ])   histKNN    ->Write();
  if (Use["HMatrix"      ])   histHm     ->Write();
  if (Use["Fisher"       ])   histFi     ->Write();
  if (Use["FisherG"      ])   histFiG    ->Write();
  if (Use["BoostedFisher"])   histFiB    ->Write();
  if (Use["LD"           ])   histLD     ->Write();
  if (Use["MLP"          ])   histNn     ->Write();
  if (Use["MLPBFGS"      ])   histNnbfgs ->Write();
  if (Use["MLPBNN"       ])   histNnbnn  ->Write();
  if (Use["CFMlpANN"     ])   histNnC    ->Write();
  if (Use["TMlpANN"      ])   histNnT    ->Write();
  if (Use["BDT"          ])   histBdt    ->Write();
  if (Use["BDTD"         ])   histBdtD   ->Write();
  if (Use["BDTG"         ])   histBdtG   ->Write();
  if (Use["RuleFit"      ])   histRf     ->Write();
  if (Use["SVM_Gauss"    ])   histSVMG   ->Write();
  if (Use["SVM_Poly"     ])   histSVMP   ->Write();
  if (Use["SVM_Lin"      ])   histSVML   ->Write();
  if (Use["FDA_MT"       ])   histFDAMT  ->Write();
  if (Use["FDA_GA"       ])   histFDAGA  ->Write();
  if (Use["Category"     ])   histCat    ->Write();
  if (Use["Plugin"       ])   histPBdt   ->Write();

  // Write also error and significance histos
  if (Use["PDEFoam"]) { histPDEFoam->Write(); histPDEFoamErr->Write(); histPDEFoamSig->Write(); }

  // Write also probability hists
  if (Use["Fisher"]) { if (probHistFi != 0) probHistFi->Write(); if (rarityHistFi != 0) rarityHistFi->Write(); }
  target->Close();

  std::cout << "--- Created root file: \"TMVA_application.root\" containing the MVA output histograms" << std::endl;

  delete reader;

  std::cout << "==> ClassificationApplication is done!" << endl << std::endl;
}
