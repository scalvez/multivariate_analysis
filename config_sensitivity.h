#ifndef conf_sensitivity_h
#define conf_sensitivity_h

namespace conf_sens {

  const double year2sec = 3600 * 24 * 365.25;

  double isotope_mass_number = 82;
  double isotope_mass = 7000; //kg
  double exposure = 2.5; //years
  double tracker_volume = 15.3; //m3

  double T_2nu = 9.0e19; //years

  double eff_0nu_2e = 24269./1e5; //uBq/kg
  double eff_2nu_2e = 2063343./1e7/25 ; //uBq/kg  // only correct above 2MeV
  double eff_tl208_2e = 10875./1e7; //uBq/kg
  double eff_bi214_2e = 29148./1e7; //uBq/kg
  double eff_radon_2e = 729./1e7; //150uBq/m3

  double tl208_activity = 2. /1e6; //Bq/kg
  double bi214_activity = 10. /1e6; //Bq/kg
  double radon_activity = 150 /1e6; //Bq/m3

  double k_sens = log(2) * 6.022e23  * isotope_mass * exposure / isotope_mass_number;

}

#endif
