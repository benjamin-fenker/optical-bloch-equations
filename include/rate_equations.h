// Authors: Benjamin Fenker 2013
// Copyright 2013 Benjamin Fenker

#ifndef INCLUDE_RATE_EQUATIONS_H_
#define INCLUDE_RATE_EQUATIONS_H_

#include <vector>
#include "./optical_pumping_method.h"
#include "./optical_pumping_data_structures.h"
#include "./eigenvector_helper.h"
// Rate equations will be based off of the formalism of
// Nafcha PRA 52 (4) 1995 Section IV.  The defining equations are 12, 13, 15.
// These are recast a 3.36 and 3.39 of DM's thesis.  DM uses \omega for his
// Lorentzian factor, presumably meaning angular frequency.  However, I beleive
// that is incorrect and will use instead frequency with E = hbar*nu.  I REALLY
// NEED TO CHECK THIS FOR SURE!!!

using std::vector;
class Rate_Equations: public OpticalPumping_Method {
 public:
  Rate_Equations();
  Rate_Equations(Eigenvector_Helper set_eigen, Laser_data laser_fe,
                 Laser_data laser_ge);
  ~Rate_Equations();
  void setup_transition_rates(double linewidth);
  double set_transition_rate(double laser_power, double sat_intensity,
                             double atom_lw, double laser_lw,
                             double atom_freq, double laser_freq);
  void update_population(double dt);
  void switch_off_laser(int las);
  static int update_population_gsl(double t, const double y[],
                               double dydt[], void *params);
  static int jacobian(double t, const double y[], double *dfdy, double dfdt[],
               void *params);
  void reset_dPop();
  vector<vector<vector<double> > > transition_rate_eg;  // e,g,q are the indexes
  vector<vector<vector<double> > > transition_rate_ef;  // e,g,q are the indexes
  vector<double> dPop_g;
  vector<double> dPop_f;
  vector<double> dPop_e;
};

#endif  // INCLUDE_RATE_EQUATIONS_H_
