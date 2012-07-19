// Copyright 2012 Benjamin Fenker


#ifndef NEWOBE_INCLUDE_OPTICAL_PUMPING_METHOD_H_
#define NEWOBE_INCLUDE_OPTICAL_PUMPING_METHOD_H_
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include "./eigenvector_helper.h"
#include "./optical_pumping_data_structures.h"

using std::vector;
// IMPORTANT:  The methods of this class are all blind to units, so you must
// pass parameters that all match!  For example MHz and ns^-1 conflict as they
// both are time^-1

class OpticalPumping_Method {
 public:
  OpticalPumping_Method();
  OpticalPumping_Method(atom_data atom, magnetic_field_data field,
                        Laser_data laser_fe, Laser_data laser_ge);
  virtual ~OpticalPumping_Method();
  void setup_quantum_numbers(int I2, int J2);
  void setup_quantum_numbers(atom_data atom);
  void setup_frequencies_excited(int I2, int Je2, double excitation,
                                 double hyperfine_const, double g_I,
                                 double mu_B, double B_z);
  void setup_frequencies_excited(atom_data atom, magnetic_field_data field);
  void setup_frequencies_ground(int I2, double hyperfine_const, double g_I,
                                double mu_B, double B_z);
  void setup_frequencies_ground(atom_data atom, magnetic_field_data field);
  double set_frequency(double excitation, int I2, int J2, int F2, int Mf2,
                       int L2, double hyperfine_const, double mu_B, double g_I,
                       double B_z);
  void setup_eg_coupling(int I2, int Je2);
  void setup_eg_coupling(atom_data atom);
  void setup_ef_coupling(int I2, int Je2);
  void setup_ef_coupling(atom_data atom);
  double set_coupling(int I2, int Jg2, int Fg2, int Mg2,
                      int Je2, int Fe2, int Me2, int q);
  void setup_pop_uniform_ground();

  double get_total_population();
  double get_polarization();
  double get_alignment();
  bool is_hermitian();


  // Virtual methods for subclasses to override
  // The =0 defines it as a pure virtual function without an implementation.
  // Any subclass MUST fully implement this function
  virtual void update_population(double dt) = 0;
  //  virtual int update_population_gsl(double t, const double y[],
  //                                double f[], void *params) = 0;
  void print_couplings(FILE * des);
  void print_density_matrix(FILE *  des);
  void print_data(FILE *des, double time);
  int numEStates, numFStates, numGStates;
  double tau;  // Natural lifetime of the excited state

  // Sturcts to hold laser polarization, intensity, frequency
  Laser_data laser_ge;
  Laser_data laser_fe;

  // Arrays to hold quantum numbers with the standard indexing
  // F = I + J ; J = L + S
  vector<int> Fe2_Vector;
  vector<int> MFe2_Vector;
  vector<int> Ff2_Vector;
  vector<int> MFf2_Vector;
  vector<int> Fg2_Vector;
  vector<int> MFg2_Vector;

  // Arrays to hold the frequency of each state.  NOT ANGULAR FREQUENCY
  // The energy difference between two states will be given by:
  // \deltaE = hbar(nu_1 - nu_2) where hbar is planck's constant
  vector<double> nu_E;
  vector<double> nu_F;
  vector<double> nu_G;

  // Arrays to hold the coupling of each state to one another.
  // See de Clerq J. Physique 45 (1984) 239-247 Eq I.6 or
  // Tremblay PRA 41 9 1990 Eq 17
  // This differs from Metcalf 1999 by a constant factor only.  This constant
  // factor is just what Metcalf needs to have his expression (4.33) normalized
  // correctly.  de Clerq equation I.2 appears different, but that could be
  // in evaluating the reduced matrix element.  Perhaps similarly with Metcalf
  vector<vector<vector<double> > > a_eg;
  vector<vector<vector<double> > > a_ef;

  /* Density matrix is broken up into siz parts to match T&J.  
     1st part is rho_ee which will hold the excited state population on the
     diagonals and excited state hyperfine and Zeeman coherences otherwise
     
     2nd part is rho_gg which is similar to rhoee except it holds the 
     F = I-1/2 ground states
     
     3rd part is rho_ff.  Same as rhogg except the F = I + 1/2 ground states
     
     4th part is delta_rho_ef which holds e-f optical coherences
     5th part is delta_rho_eg which holds e-g optical coherences
     6th part is delta_rho_fg which holds ground state hyperfine coherences
  */
  vector<vector<gsl_complex> > rho_ee;
  vector<vector<gsl_complex> > rho_ff;
  vector<vector<gsl_complex> > rho_gg;
  vector<vector<gsl_complex> > delta_ef;
  vector<vector<gsl_complex> > delta_eg;
  vector<vector<gsl_complex> > delta_fg;
};
#endif  // NEWOBE_INCLUDE_OPTICAL_PUMPING_METHOD_H_
