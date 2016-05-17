#include <iostream>
#include "channel_selection.h"
#include "pseudo_generator.h"
#include <TTree.h>
#include <TApplication.h>
#include "sensitivity_measurements.h"

int main() {

  //Required to prevent crash : force the loading of all required library like TTree
  TApplication *myapp=new TApplication("myapp",0,0);

  //change to map where <isotope,activity>
  // std::vector<std::string> isotopes = {"tl208", "bi214"};
  // std::vector<double> activity = {100e-6, 100e-6};
  std::map < std::string, double > isotope_activity;
  isotope_activity.insert(std::pair<std::string,double>("tl208",100e-6));
  isotope_activity.insert(std::pair<std::string,double>("bi214",100e-6));

  //maybe later moved to analysis_config.h
  std::vector <TString> quantities;
  quantities.push_back("1e1g_electron_gamma_energy_sum");
  quantities.push_back("1e2g_electron_gammas_energy_sum");

  for (auto i = isotope_activity.begin(); i != isotope_activity.end(); ++i) {
    const std::string & a_isotope = i->first;
    const double & a_activity = i->second;
    channel_selection(a_isotope, quantities, quantity_efficiency);
  }

  for (auto i = isotope_activity.begin(); i != isotope_activity.end(); ++i) {
    const std::string & a_isotope = i->first;
    const double & a_activity = i->second;
    pseudo_generator(a_isotope,quantities,a_activity,quantity_efficiency);
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
