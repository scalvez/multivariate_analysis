#ifndef ANALYSIS_CONFIG_H
#define ANALYSIS_CONFIG_H 1

#include "TCut.h"
#include "TString.h"

// const double bi214_channel_1e1a_efficiency = 10733./2500000;

// const double tl208_channel_1e1g_efficiency = 0.139637;
// const double bi214_channel_1e1g_efficiency = 0.123019;

// 2e_int channel
const double internal_probability_min = 0;
const double internal_probability_max = 0;
const double vertices_probability_min = 0;
const double vertices_probability_max = 0;

// TCut prob_int_cut = "2e_electrons_internal_probability > 0.04";
// TCut vertices_proba_cut = "2e_electrons_vertices_probability > 0.04";
// TCut channel_2e_int_cut = prob_int_cut&&vertices_proba_cut;

void get_histogram_options(TString quantity, int & nbins, double & xmin, double & xmax);

int get_isotope_mc_size(TString isotope);

#endif
