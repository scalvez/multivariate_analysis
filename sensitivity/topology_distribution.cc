#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <vector>

void topology_distribution()
{
  TString input_file = "./trees/2nu_tree.root";
  TString output_file = "./topo_distrib/2nu_topo_distrib.root";
  // TString input_file = "./trees/tl208_tree.root";
  // TString output_file = "./topo_distrib/tl208_topo_distrib.root";

  double mc_size = 1e8;

  TFile *f = TFile::Open(input_file);
  TTree *tree = (TTree*)f->Get("snemodata");
  TFile *f_output= new TFile(output_file, "RECREATE");

  std::vector<TString> quantities_pdf;
  // quantities_pdf.push_back("1e_electron_energy");
  // quantities_pdf.push_back("1e1a_electron_energy");
  // quantities_pdf.push_back("2e_electrons_energy_sum");
  // quantities_pdf.push_back("1e1p_electron_energy");
  // quantities_pdf.push_back("2p_positrons_energy_sum");
  // quantities_pdf.push_back("1e1g_electron_energy");
  // quantities_pdf.push_back("1e2g_electron_energy");
  // quantities_pdf.push_back("1e3g_electron_energy");
  quantities_pdf.push_back("2e1g_electrons_gammas_energy_sum");
  // quantities_pdf.push_back("2e2g_electrons_gammas_energy_sum");
  // quantities_pdf.push_back("2e3g_electrons_gammas_energy_sum");

  double other_topologies = mc_size;

  TH1F *h_td = new TH1F("h_topology_dsitribution","",12,1,13);
  for(unsigned int j=0; j<quantities_pdf.size(); ++j) {

    int nbins=100;
    double xmin=0, xmax=5;

    TString qty = quantities_pdf.at(j);

    TH1F* h = new TH1F(qty,qty,nbins,xmin,xmax);

    /*
    if(qty.Contains("2e1g")) {
      std::cout << "debug 2e1g topology " << std::endl;
      TCut cut = "2e1g_electrons_gammas_energy_sum > 0.01";
      tree->Project(qty,qty,cut);
    }
    else if(qty.Contains("2e2g")) {
      TCut cut = "2e2g_electrons_gammas_energy_sum != 0";
      tree->Project(qty,qty,cut);
    }
    else if(qty.Contains("2e3g")) {
      TCut cut = "2e3g_electrons_gammas_energy_sum != 0";
      tree->Project(qty,qty,cut);
    }
    else
      tree->Project(qty,qty);
    */

    TCut cut = "2e1g_electron_minimal_energy != 0 || 2e1g_gamma_energy != 0";
    tree->Project(qty,qty,cut);
    // TCut cut = "2e1g_electron_minimal_energy != 0";
    // tree->Project(qty,qty,cut);
    // tree->Project(qty,qty);

    std::cout << h->Integral(1,100) << std::endl;
    /*
    // h->ClearUnderflowAndOverflow(); Unavailable with this version of ROOT
    h->SetBinContent(0,0);
    h->SetBinContent(nbins+1,0);
    other_topologies -= h->Integral(1,100);
    h_td->Fill(j+1,h->Integral(1,100)/mc_size);
    TString topo(qty(0, qty.First("_")));

    // h_td->GetXaxis()->SetBinLabel(j+1, Form("(%d) (%d)",sta1));
    double test = h->Integral(1,100)/mc_size;

    TString label = topo + "\n" + Form(" %d",double(h->Integral(1,100)/mc_size));
    std::cout << "label " << label << std::endl;
    std::cout << "  eff " << test << std::endl;
    h_td->GetXaxis()->SetBinLabel(j+1, label);
    // std::cout << topo << "  " << h->Integral(1,100)/mc_size << std::endl;
    h_td->Draw();
    */
  }

  if(other_topologies > 0) {
    h_td->Fill(12,other_topologies/mc_size);
    h_td->GetXaxis()->SetBinLabel(12, "other \n");
    std::cout << "other" << "  " << other_topologies/mc_size << std::endl;
  }

  // h_td->GetXaxis()->LabelsOption("v");

  h_td->GetXaxis()->SetTitle("Topology");
  h_td->GetXaxis()->SetTitleFont(62);
  h_td->GetXaxis()->SetTitleSize(0.05);
  h_td->GetXaxis()->SetTitleOffset(0.88);
  h_td->GetXaxis()->SetRangeUser(0.,14.);
  h_td->GetYaxis()->SetTitle("Efficiency");
  h_td->GetYaxis()->SetTitleFont(62);
  h_td->GetYaxis()->SetTitleSize(0.05);
  h_td->GetYaxis()->SetTitleOffset(0.9);
  h_td->GetYaxis()->SetTitleOffset(0.9);
  h_td->GetYaxis()->SetRangeUser(0.00001,1.);
  h_td->SetMarkerStyle(20);
  h_td->SetMarkerSize(1);
  h_td->SetMarkerColor(kGreen+1);
  h_td->SetLineColor(kGreen+1);
  h_td->SetFillColor(kGreen+1);
  h_td->SetDrawOption("A");

  h_td->Write();

  return;
}
