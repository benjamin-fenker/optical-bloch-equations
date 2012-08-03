// Copyright Benjamin Fenker 2012

#ifndef INCLUDE_EIGENVECTOR_HELPER_H_
#define INCLUDE_EIGENVECTOR_HELPER_H_

#include "./optical_pumping_data_structures.h"

class Eigenvector_Helper {
 public:
  void diagH(int I2, int L, int J2, double mu_B, double B_z, double B_x,
             double g_I, double Aj, double *arrayToFill);
  void diagH(atom_data atom, magnetic_field_data field, int L, int J2,
             double Aj, double *arrayToFill);
  double calc_gj(int J2, int L2, int S2);
  double calc_gf(int F2, int J2, int I2, int L2, int S2, double g_I);
};
#endif  // INCLUDE_EIGENVECTOR_HELPER_H_
