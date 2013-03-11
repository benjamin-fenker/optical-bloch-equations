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
  void calculate_derivs(DM_container *status);
  /* static int update_population_gsl(double t, const double y[], double f[], */
  /*                                  void *params); */

  void setup_dipole_moments(double gamma);
  double set_dipole_moment(double gamma, double omega_ex, double omega_gr);
  void print_rabi_frequencies(FILE * des);
  /* void reset_dPop(); */


  void integrate_ee(DM_container *in);
  void integrate_gg(DM_container *in);
  void integrate_ff(DM_container *in);
  void integrate_eg(DM_container *in);
  void integrate_ef(DM_container *in);
  void integrate_fg(DM_container *in);

  void change_magnetic_field(double newfield);
  vector<vector<double> > dipole_moment_eg;
  vector<vector<double> > dipole_moment_ef;

  /* vector<vector<gsl_complex> > drho_ee; */
  /* vector<vector<gsl_complex> > drho_ff; */
  /* vector<vector<gsl_complex> > drho_gg; */
  /* vector<vector<gsl_complex> > ddelta_ef; */
  /* vector<vector<gsl_complex> > ddelta_eg; */
  /* vector<vector<gsl_complex> > ddelta_fg; */

  bool es_Zeeman, gs_Zeeman, es_hyperfine, gs_hyperfine;
};

#endif  // INCLUDE_DENSITY_MATRIX_H_
