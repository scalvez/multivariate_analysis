#include "analysis_config.h"

// TCut get_channel_cut (TSring channel) {

//   switch(channel)
//     {
//     case "2e_int" :
//       return channel_2e_int_cut;
//       break;
//     case "1e1a" :
//       return channel_1e1a_cut;
//       break;
//     }
//   return "";
// }

// void get_histogram_options(TString quantity, int & nbins, double & xmin, double & xmax) {

//   // hardcoded for now, see if it is parametrized in the configuration file
//   if(quantity.Contains("energy")) {
//     nbins = 100;
//     xmin = 0;
//     xmax = 5;
//   }
//   else if(quantity.Contains("probability")){
//     nbins = 100;
//     xmin = 0;
//     xmax = 1;
//   }
//   else if(quantity.Contains("angle")) {
//     nbins = 100;
//     xmin = -1;
//     xmax = 1;
//   }
//   else if (quantity.Contains("track_length")) {
//     nbins = 100;
//     xmin = 0;
//     xmax = 500;
//   }
//   else {
//     //Also maybe the alpha delayed time
//     //for now, dummy values
//     nbins = 100;
//     xmin = 0;
//     xmax = 1000;
//   }

//   return;
// }
