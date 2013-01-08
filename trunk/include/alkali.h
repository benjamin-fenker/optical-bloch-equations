// Authors: Benjamin Fenker 2013
// Copyright 2012 Benjamin Fenker

#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include "./optical_pumping_data_structures.h"
#include "./rate_equations.h"

#ifndef INCLUDE_ALKALI_H_
#define INCLUDE_ALKALI_H_

using std::string;

class Alkali {
 public:
  int lookupParameters(string isotope, int Je2, int* I2, double* Aj_g,
                       double* Aj_e, double* g_I,
                       double* omega_Di, double *tau);

  int getNumberOfExcitedStates(int I2, int Je2);
  int getNumberOfGroundStates_f(int I2);
  int getNumberOfGroundStates_g(int I2);

  double getLaserFrequency(atom_data atom, magnetic_field_data field, int Fg2,
                           int Mfg2, int Fe2, int Mfe2, double detune);
  double getGamma(double tau, double laser_power);
};

#endif  // INCLUDE_ALKALI_H_
