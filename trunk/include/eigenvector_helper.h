// Copyright Benjamin Fenker 2012

#ifndef INCLUDE_EIGENVECTOR_HELPER_H_
#define INCLUDE_EIGENVECTOR_HELPER_H_

#include <vector>
#include "./optical_pumping_data_structures.h"

using std::vector;

class Eigenvector_Helper {
 public:
  Eigenvector_Helper() {}
  Eigenvector_Helper(atom_data set_atom, magnetic_field_data set_field);

  /*
  void diagH(int I2, int L, int J2, double mu_B, double B_z, double B_x,
             double g_I, double Aj, double *arrayToFill);
  void diagH(atom_data atom, magnetic_field_data field, int L, int J2,
             double Aj, double *arrayToFill);
  */
  vector<vector<double> > diagH(int L);

  //  void get_nuclear_spin_vec(int numBasisStates, int I2, int *Iz2);
  vector<int> get_nuclear_spin(int numBasisStates);
  vector<int> get_total_atomic_spin(int numBasisStates, int J2);
  //  void get_total_atomic_spin_vec(int numBasisStates, int I2, int J2, int *Jz2);
  double calc_gj(int J2, int L2, int S2);
  double calc_gf(int F2, int J2, int I2, int L2, int S2, double g_I);

  atom_data atom;
  magnetic_field_data field;

  // Leave room for odd indices to hold imaginary components
  vector<int> nuclear_spin_ground;
  vector<int> nuclear_spin_excited;
  vector<int> total_atomic_spin_ground;
  vector<int> total_atomic_spin_excited;
  vector<vector<double> > IzJz_decomp_ground;
  vector<vector<double> > IzJz_decomp_excited;
};
#endif  // INCLUDE_EIGENVECTOR_HELPER_H_
