#include "analysis_config.h"

TCut get_channel_cut (TSring channel) {

  switch(channel)
    {
    case "2e_int" :
      return channel_2e_int_cut;
      break;
    case "1e1a" :
      return channel_1e1a_cut;
      break;
    }
  return "";
}
