// Authors: Benjamin Fenker 2013
// Copyright 2013 Benjamin Fenker


#ifndef INCLUDE_OPTICAL_PUMPING_METHOD_H_
#define INCLUDE_OPTICAL_PUMPING_METHOD_H_
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
  OpticalPumping_Method(Eigenvector_Helper set_eigen, Laser_data laser_fe,
                        Laser_data laser_ge, double tilt);
  virtual ~OpticalPumping_Method();
  void setup_quantum_numbers(int I2, int J2);
  void setup_quantum_numbers(atom_data atom);
  void setup_gFactors(atom_data atom);
  void setup_frequencies_excited(int I2, int Je2, double excitation,
                                 double hyperfine_const,
                                 double B_z);
  void setup_frequencies_excited(atom_data atom, magnetic_field_data field);
  void setup_frequencies_ground(int I2, double hyperfine_const,
                                double B_z);
  void setup_frequencies_ground(atom_data atom, magnetic_field_data field);
  /* const guarantees that this won't change the state */
  double set_frequency(double excitation, int I2, int J2, int F2, int Mf2,
                       int L2, double hyperfine_const, double g_I,
                       double B_z);
  double set_frequency(double excitation, int I2, int J2, int F2, int Mf2,
                       double hyperfine_const, double B_z, double _gf);
  void setup_eg_coupling(int I2, int Je2);
  void setup_eg_coupling(atom_data atom);
  void setup_ef_coupling(int I2, int Je2);
  void setup_ef_coupling(atom_data atom);
  double set_coupling(int I2, int Jg2, int Fg2, int Mg2,
                      int Je2, int Fe2, int Me2, int q);
  void setup_raising();
  void setup_lowering();
  double set_raising(int F2, int M2);
  double set_lowering(int F2, int M2);
  void setup_pop_uniform_ground();
  void setup_pop_withTilt(double tilt);
  virtual void switch_off_laser(int las); /* las = 1: g-->e; las = 2: f-->e */
  virtual void change_magnetic_field(double newField);
  // Struct to hold atom, laser and eignenvector information
  Eigenvector_Helper eigen;

  double get_total_population();
  double get_polarization();
  double get_alignment();
  double get_excited_state_total();
  bool is_hermitian();

  op_data_for_gsl data;                 /* Not really for GSL anymore */

  //  store the entire density matrix
  int totalTerms;

  void reset_dPop();
  void apply_dPop(double dt);

  // The =0 defines it as a pure virtual function without an implementation.
  // Any subclass MUST fully implement this function
  virtual void calculate_derivs(DM_container *status) = 0;

  void update_population_euler(double dt);
  void update_population_RK4(double dt);

  void apply_transverse_field(DM_container *in);
  void apply_transverse_field_ee(DM_container *in);
  void apply_transverse_field_gg(DM_container *in);
  void apply_transverse_field_ff(DM_container *in);
  void apply_transverse_field_eg(DM_container *in);
  void apply_transverse_field_ef(DM_container *in);
  void apply_transverse_field_fg(DM_container *in);

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
  // \deltaE = h(nu_1 - nu_2) where h is planck's constant
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
  /* Gyromagetic factor for each state */
  vector<double> gFactor_E;     /* Different g-Factor for each F-manifold    */
  double gFactor_F;             /* All have F = I + 1/2 so g-factor the same */
  double gFactor_G;             /* All have F = I - 1/2 so g-factor the same */

  /* Raising and lowering eigenvalues for each state
     (see Renzoni 2001, below 2c.  Store them here so don't have to calculate
     them at every time step */
  vector<double> cPlus_E;
  vector<double> cPlus_F;
  vector<double> cPlus_G;

  vector<double> cMins_E;
  vector<double> cMins_F;
  vector<double> cMins_G;

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


  /* vector<vector<gsl_complex> > rho_ee; */
  /* vector<vector<gsl_complex> > rho_ff; */
  /* vector<vector<gsl_complex> > rho_gg; */
  /* vector<vector<gsl_complex> > delta_ef; */
  /* vector<vector<gsl_complex> > delta_eg; */
  /* vector<vector<gsl_complex> > delta_fg; */
  DM_container *dm_status, *dm_derivs;
  /* Matrices to hold the ground and excited state decomposition. */
  vector<vector<double> > IzJz_ground;
  vector<vector<double> > IzJz_excited;

  bool es_Zeeman, gs_Zeeman, es_hyperfine, gs_hyperfine;
};
#endif  // INCLUDE_OPTICAL_PUMPING_METHOD_H_
