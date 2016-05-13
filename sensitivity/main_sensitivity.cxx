#include <iostream>
#include "channel_selection.h"

int main() {

  std::vector<std::string> isotopes = {"tl208", "bi214"};

  std::vector <TString> isotope_mc_files;
  std::vector <TString> isotope_pdf_files;
  // std::map <TString,double> quantities_pdf;
  std::vector <TString> quantities_pdf;

  quantities_pdf.push_back("1e1g_electron_gamma_energy_sum");

  for (auto i = isotopes.begin(); i != isotopes.end(); ++i) {
    const std::string & a_isotope = *i;

    isotope_mc_files.push_back("../" + a_isotope + "_tree.root");
  // for(unsigned int i = 0; i < quantities_pdf.size(); ++i)
  //   quantity_efficiency[i];
    isotope_pdf_files.push_back("../" + a_isotope + "_pdf.root");
  }

  channel_selection(isotope_mc_files, isotope_pdf_files, quantities_pdf);

  TString pseudo_file = "pseudo.root";

  // std::cout << "vector size is " << isotope_quantity_efficiency.size() << std::endl;
  // std::map < TString, double > test = isotope_quantity_efficiency.at(0);
  // std::cout << "map size is " << test.size() << std::endl;

  // std::map<TString, double>::iterator iter = test.begin();
  // std::cout << (*iter).first << std::endl;

  // std::vector<TString>::iterator iter = isotope_mc_files.begin();

  // std::cout << *iter << std::endl;

  // std::cout << " first  " << (*it).first << "    second  " << (*it).second << std::endl;
  // std::cout << test["1e1g_electron_gamma_energy_sum"] << std::endl;
}
