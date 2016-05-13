#ifndef CHANNEL_SELECTION_H
#define CHANNEL_SELECTION_H 1

#include <TString.h>
#include <vector>
#include <map>

void get_histogram_options(TString quantity, int & nbins, double & xmin, double & xmax);

double get_isotope_mc_size(TString isotope);

void channel_selection(std::vector <TString> input_files);// , std::vector<TString> output_files, std::vector<TString> quantities_pdf, std::vector < std::map < TString,double > > & isotope_qty_eff, bool normalize = true);

//void channel_selection();
#endif
