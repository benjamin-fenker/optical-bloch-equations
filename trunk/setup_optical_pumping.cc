// Copyright 2012 Benjamin Fenker

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <iomanip>
#include "include/optical_pumping.h"
#include "include/units.h"

using std::string;

bool op_verbose = false;

int main(int argc, char* argv[]) {
  printf("\n");

  double tmax = 1.0 * _ns;  // ns
  double dt = 1.0 * _ns;  // ns

  string method = "R";  // O for OBE and R fro Rate Equations

  bool zCoherences = false;
  bool hfCoherences_ex = false;
  bool hfCoherences_gr = false;

  string isotope = "37K";
  int Je2 = 1;  // D1(1) or D2(3) line
  double laser_fe_I = 800.0 * (_uW/_cm2);  // uW/cm^2
  double laser_ge_I = 800.0 * (_uW/_cm2);  // uW/cm^2

  // (I think) that John is reportings I_+ + I_- / I_+ - I_-, which is
  // equivalent to s3/s0 in the notation of Jackson 3rd edition
  double laser_fe_s3_over_s0 = -0.995; // 
  double laser_ge_s3_over_s0 = -0.995; // 

  double laser_fe_detune = -4.5 * _MHz;  // MHz
  double laser_ge_detune = -4.5 * _MHz;  // MHz

  double laser_fe_linewidth = 0.2 *_MHz;  // MHz (FWHM)
  double laser_ge_linewidth = 0.2 *_MHz;  // MHz (FWHM)

  double B_z = 2.0 * _G;  // G

  // Also will accept command line input
  if (argc > 1) {
    if (strcmp(argv[1], "-h") == 0) {
      printf("First paramter is isotope (37K, 38K, 47K) [%s]\n",
             isotope.c_str());
      printf("Second parameter is Je2 (1 for D1 line, 3 for D2 line [%d]\n",
             Je2);
      printf("Third parameter is tmax in ns [%3.1G]\n", tmax/_ns);
      printf("Fourth paramter is dt in ns [%3.1G]\n", dt/_ns);
      printf("Fifth parameter is B_ext in G [%3.1G]\n", B_z/_G);
      printf("Sixth parameter is detune (both lasers) in MHz [%3.1G]\n",
             laser_fe_detune/_MHz);
      printf("Seventh paramter is linewidth (both lasers) in MHz [%3.1G]\n",
             laser_fe_linewidth/_MHz);
      printf("\n\n");
      return 0;
    } else {
      isotope = argv[1];
      if (argc>2) {
        Je2 = atoi(argv[2]);
        if (argc > 3) {
          tmax = atof(argv[3])*_ns;
          if (argc > 4) {
            dt = atof(argv[4])*_ns;
            if (argc > 5) {
              B_z = atof(argv[5])*_G;
              if (argc > 6) {
                laser_fe_detune = atof(argv[6]) *_MHz;
                laser_ge_detune = atof(argv[6]) *_MHz;
                if (argc > 7) {
                  laser_fe_linewidth = atof(argv[7]) *_MHz;
                  laser_ge_linewidth = atof(argv[7]) *_MHz;
                }
              }
            }
          }
        }
      }
    }
  }

  OpticalPumping pumper;
  int status = pumper.pump(isotope, method, tmax, dt, zCoherences,
                           hfCoherences_ex, hfCoherences_gr, Je2, laser_fe_I,
                           laser_ge_I, laser_fe_detune, laser_ge_detune,
                           laser_fe_linewidth, laser_ge_linewidth,
                           laser_fe_s3_over_s0, laser_ge_s3_over_s0, B_z);
  printf("\nCompleted with status = %d\n\n", status);
  return status;
}
