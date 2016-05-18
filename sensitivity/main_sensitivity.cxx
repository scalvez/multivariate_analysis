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

  for (auto i = isotope_activity.begin(); i != isotope_activity.end(); ++i) {
    const std::string & a_isotope = i->first;
    const double & a_activity = i->second;
    if (generate_pdf)
      channel_selection(a_isotope, quantities, quantity_efficiency);
    if (generate_pseudo)
      pseudo_generator(a_isotope,quantities,a_activity, quantity_efficiency);
  }

  gSystem->Exec("hadd -f ../pseudo.root ../*_pseudo.root");

  // for (auto i = isotope_activity.begin(); i != isotope_activity.end(); ++i) {
  //   const std::string & a_isotope = i->first;
  //   const double & a_activity = i->second;
  //   if (generate_pseudo)
  //     pseudo_generator(a_isotope,quantities,a_activity,quantity_efficiency);
  //   const std::map < std::string, double > isotope_activity = {
  //     {"tl208",100e-6},
  //     {"bi214",100e-6}
  //   };
  // }

  std::map < std::string, std::vector<double> > activity_measurement;
  multi_fit(activity_measurement);

  std::cout << " - Reconstructed activities - " << std::endl;
  for (auto i = activity_measurement.begin(); i != activity_measurement.end(); ++i) {
    std::string isotope = i->first;
    double activity = i->second.at(0);
    std::cout << " Isotope : " << isotope << std::endl;
    std::cout << " Activity : " << activity * 1e6  << " uBq/kg " << std::endl;
  }

  quantity_efficiency.clear();
  return 0;
}

void main_sensitivity() {
  main();
  return;
}
