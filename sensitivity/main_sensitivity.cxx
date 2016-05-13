#include <iostream>
#include "channel_selection.h"

int main() {

  // std::vector<std::string> isotopes = {"tl208", "bi214"};
  std::vector<std::string> isotopes;
  isotopes.push_back("tl208");

  //maybe later moved to analysis_config.h
  std::vector <TString> quantities_pdf;
  quantities_pdf.push_back("1e1g_electron_gamma_energy_sum");
  quantities_pdf.push_back("1e2g_electron_gammas_energy_sum");

  for (auto i = isotopes.begin(); i != isotopes.end(); ++i) {
    std::cout << " " << *i << std::endl;
    const std::string & a_isotope = *i;
    channel_selection(a_isotope, quantities_pdf);
  }

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

void main_sensitivity() {
  main();
  return;
}
