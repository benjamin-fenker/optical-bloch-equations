// Copyright 2012 Benjamin Fenker

#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include "./optical_pumping_data_structures.h"

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

  double getLaserFrequency(double omega_Di, int Fg2, int I2, int Je2, int Fe2,
                           double Aj_g, double Aj_e, double detune);
  double getLaserFrequency(atom_data atom, int Fg2, int Fe2, double detune);
  double getGamma(double tau, double laser_power);
};

#endif  // INCLUDE_ALKALI_H_
