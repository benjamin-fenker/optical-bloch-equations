// Authors: Benjamin Fenker 2013
// Copyright 2012 Benjamin Fenker

#ifndef INCLUDE_DENSITY_MATRIX_H_
#define INCLUDE_DENSITY_MATRIX_H_

#include <vector>
#include "./optical_pumping_method.h"
#include "./optical_pumping_data_structures.h"
#include "./eigenvector_helper.h"

using std::vector;

class Density_Matrix: public OpticalPumping_Method {
 public:
  Density_Matrix(Eigenvector_Helper set_eigen, Laser_data laser_fe,
                 Laser_data laser_ge, coherence_flags flags);
  //  ~Density_Matrix();
  void update_population(double dt);
  static int update_population_gsl(double t, const double y[], double f[],
                                   void *params);

  void setup_dipole_moments(double gamma);
  double set_dipole_moment(double gamma, double omega_ex, double omega_gr);
  void print_rabi_frequencies(FILE * des);
  void reset_dPop();
  void apply_dPop(double dt);

  void integrate_ee();
  void integrate_gg();
  void integrate_ff();
  void integrate_eg();
  void integrate_ef();
  void integrate_fg();

  void change_magnetic_field(double newfield);
  vector<vector<double> > dipole_moment_eg;
  vector<vector<double> > dipole_moment_ef;

  vector<vector<gsl_complex> > drho_ee;
  vector<vector<gsl_complex> > drho_ff;
  vector<vector<gsl_complex> > drho_gg;
  vector<vector<gsl_complex> > ddelta_ef;
  vector<vector<gsl_complex> > ddelta_eg;
  vector<vector<gsl_complex> > ddelta_fg;

  bool es_Zeeman, gs_Zeeman, es_hyperfine, gs_hyperfine;
};

#endif  // INCLUDE_DENSITY_MATRIX_H_
