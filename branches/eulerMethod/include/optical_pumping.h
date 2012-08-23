// Copyright 2012 Benjamin Fenker

#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <string>

#ifndef INCLUDE_OPTICAL_PUMPING_H_
#define INCLUDE_OPTICAL_PUMPING_H_

using std::string;

class OpticalPumping {
 public:
  int pump(string isotope, string method, double tmax, double tStep,
           bool zCoherences, bool hfCoherences_ex, bool hfCoherences_gr,
           int Je2, double laser_fe_power, double laser_ge_power,
           double laser_fe_detune, double laser_ge_detune,
           double laser_fe_linewidth, double laser_ge_linewidth,
           double laser_fe_pol[], double laser_ge_pol[], double B_z);
};
#endif  // INCLUDE_OPTICAL_PUMPING_H_
