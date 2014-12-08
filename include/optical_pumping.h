// Authors: Benjamin Fenker 2013
// Copyright 2013 Benjamin Fenker

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
           int Je2, int nominalSublevelTune2_ef, int nominalSublevelTune2_eg,
           double laser_fe_power, double laser_ge_power,
           double laser_fe_detune, double laser_ge_detune,
           double laser_fe_linewidth, double laser_ge_linewidth,
           double laser_fe_s3_over_s0, double laser_ge_s3_over_s0,
           double laser_fe_offTime, double laser_ge_offTime, double B_z,
           double B_x, string outfile, double tilt, int verbose);
  int test();
};
#endif  // INCLUDE_OPTICAL_PUMPING_H_

