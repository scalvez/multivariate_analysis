#include <iostream>
#include <TTree.h>
#include <TSystem.h>
#include <TApplication.h>
#include "channel_selection.h"
#include "pseudo_generator.h"
#include "multi_fit.h"
#include "sensitivity_measurements.h"
#include "analysis_config.h"

int main() {

  //Required to prevent crash : force the loading of all required library like TTree
  TApplication *myapp=new TApplication("myapp",0,0);

  for (auto i = isotope_activity.begin(); i != isotope_activity.end(); ++i) {
    const std::string & a_isotope = i->first;
    const double & a_activity = i->second;
    std::cout << "debug channel selection " << a_isotope << "  " << a_activity << std::endl;
    if (generate_pdf)
      channel_selection(a_isotope, quantities, quantity_efficiency);
    // channel_selection(a_isotope, quantities);
    pseudo_generator(a_isotope,quantities,a_activity, quantity_efficiency);
  }
  std::cout << "main ended " << std::endl;

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

  //   std::map < std::string, std::vector<double> > activity_measurement;
  //   multi_fit(activity_measurement);
  // }

}

void main_sensitivity() {
  main();
  return;
}
