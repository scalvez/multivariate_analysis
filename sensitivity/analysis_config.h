#ifndef ANALYSIS_CONFIG_H
#define ANALYSIS_CONFIG_H 1

#include "TCut.h"
#include "TString.h"
#include <vector>
#include <string>
#include <map>

const bool generate_pdf = true;
const bool generate_pseudo = true;
const bool poisson_pseudo = false;
const bool fit = true;

// const double bi214_channel_1e1a_efficiency = 10733./2500000;

// const double tl208_channel_1e1g_efficiency = 0.139637;
// const double bi214_channel_1e1g_efficiency = 0.123019;

// 2e_int channel beta-beta-like events
const double internal_probability_min = 0;
const double internal_probability_max = 0;
const double vertices_probability_min = 0;
const double vertices_probability_max = 0;

// TCut good_internal_probability_cut = "2e_electrons_internal_probability > 0.04";
// TCut vertices_proba_cut = "2e_electrons_vertices_probability > 0.04";
// TCut channel_2e_int_cut = prob_int_cut&&vertices_proba_cut;

const std::map < std::string, double > isotope_activity = {
  {"tl208",20e-6},
  {"bi214",100e-6}
};

const std::vector <TString> quantities = {"1e1g_electron_gamma_energy_sum", "1e2g_electron_gammas_energy_sum"};

void get_histogram_options(TString quantity, int & nbins, double & xmin, double & xmax);

int get_isotope_mc_size(TString isotope);

#endif
