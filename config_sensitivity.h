#ifndef conf_sensitivity_h
#define conf_sensitivity_h

namespace conf_sens {

  const double year2sec = 3600 * 24 * 365.25;

  const double isotope_mass_number = 82;
  const double isotope_mass = 7; //kg
  const double exposure = 2.5; //years
  const double tracker_volume = 15.3; //m3

  const double T_2nu = 9.0e19; //years

  const double eff_0nu_2e = 24269./1e5; //uBq/kg
  const double eff_2nu_2e = 2063343./1e7/25 ; //uBq/kg  // only correct above 2MeV
  const double eff_tl208_2e = 10875./1e7; //uBq/kg
  const double eff_bi214_2e = 29148./1e7; //uBq/kg
  const double eff_radon_2e = 729./1e7; //150uBq/m3

  const double eff_0nu_1e = 6831./1e5; //uBq/kg
  const double eff_2nu_1e = 1039397./1e7/25 ; //uBq/kg  // only correct above 2MeV
  const double eff_tl208_1e = 231932./1e7; //uBq/kg
  const double eff_bi214_1e = 1230205./1e7; //uBq/kg
  const double eff_radon_1e = 104205./1e7; //150uBq/m3

  const double tl208_activity = 2. /1e6; //Bq/kg
  // const double tl208_activity = 23.6 /1e6; //Bq/kg
  const double bi214_activity = 10. /1e6; //Bq/kg
  // const double bi214_activity = 338. /1e6; //Bq/kg
  const double radon_activity = 150 /1e6; //Bq/m3

  const double k_sens = log(2) * 6.022e23  * isotope_mass * 1000 * exposure / isotope_mass_number;

}

#endif
