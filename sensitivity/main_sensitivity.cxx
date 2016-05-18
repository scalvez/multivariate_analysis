#include <iostream>
#include <string>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>
#include "channel_selection.h"
#include "pseudo_generator.h"
#include "multi_fit.h"
#include "sensitivity_measurements.h"
#include "analysis_config.h"

int main() {

  std::cout << " Starting the sensitivity measurement..." << std::endl;

  std::cout << " Channel selection..." << std::endl;
  if (generate_pdf) {
    for (auto i = isotope_activity.begin(); i != isotope_activity.end(); ++i) {
      const std::string & a_isotope = i->first;
      std::cout << "      Isotope : " << a_isotope << std::endl;
      channel_selection(a_isotope, quantities);
    }
  }

  std::cout << " Channel selection operated" << std::endl;

  std::cout << " Pseudo experiments generation " << std::endl;
  unsigned int number_of_pseudo_experiments = 1;
  for(unsigned int n_pseudo = 0; n_pseudo < number_of_pseudo_experiments; ++n_pseudo) {
    std::map < std::string, std::vector<double> > activity_measurement;
    for (auto i = isotope_activity.begin(); i != isotope_activity.end(); ++i) {
      const std::string & a_isotope = i->first;
      const double & a_activity = i->second;
      if (generate_pseudo)
        pseudo_generator(a_isotope,quantities,a_activity);
    }
    gSystem->Exec("hadd -f ../pseudo.root ../*_pseudo.root");

    multi_fit(activity_measurement);
    std::cout << " - Reconstructed activities - " << std::endl;
    for (auto i = activity_measurement.begin(); i != activity_measurement.end(); ++i) {
      std::string isotope = i->first;
      double activity = i->second.at(0);
      std::cout << " Isotope : " << isotope << std::endl;
      std::cout << " Activity : " << activity * 1e6  << " uBq/kg " << std::endl;
    }
  }

  quantity_efficiency.clear();
  return 0;
}

void main_sensitivity() {
  main();
  return;
}
