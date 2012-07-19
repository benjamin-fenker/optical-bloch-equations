// Copyright 2012 Benjamin Fenker

#ifndef NEWOBE_INCLUDE_DENSITY_MATRIX_H_
#define NEWOBE_INCLUDE_DENSITY_MATRIX_H_

#include <vector>
#include "./optical_pumping_method.h"
#include "./optical_pumping_data_structures.h"

using std::vector;

class Density_Matrix: public OpticalPumping_Method {
 public:
  Density_Matrix(atom_data atom, magnetic_field_data field, Laser_data laser_fe,
                 Laser_data laser_ge, coherence_flags flags);
  void update_population(double dt);
  // int update_population_gsl(double t, const double y[],
  //                         double f[], void *params);

  void setup_dipole_moments(double gamma);
  double set_dipole_moment(double gamma, double omega_ex, double omega_gr);
  void print_rabi_frequencies(FILE * des);
  void reset_dPop();

  vector<vector<double> > dipole_moment_eg;
  vector<vector<double> > dipole_moment_ef;

  vector<vector<gsl_complex> > drho_ee;
  vector<vector<gsl_complex> > drho_ff;
  vector<vector<gsl_complex> > drho_gg;
  vector<vector<gsl_complex> > ddelta_ef;
  vector<vector<gsl_complex> > ddelta_eg;
  vector<vector<gsl_complex> > ddelta_fg;

  bool es_Zeeman, gs_Zeeman, es_hyperfine;
};

#endif  // NEWOBE_INCLUDE_DENSITY_MATRIX_H_
