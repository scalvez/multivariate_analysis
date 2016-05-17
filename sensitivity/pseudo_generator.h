#ifndef PSEUDO_GENERATOR_H
#define PSEUDO_GENERATOR_H 1

#include <TString.h>
#include <vector>
#include <map>

// void pseudo_generator(TString isotope, std::vector<TString> quantities, double activity, std::map < TString , double > & quantity_efficiency);
void pseudo_generator(TString isotope, std::vector<TString> quantities, double activity);

#endif
