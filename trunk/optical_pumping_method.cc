// Authors: Benjamin Fenker 2013
// Copyright 2012 Benjamin Fenker

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_coupling.h>
#include "include/optical_pumping_method.h"
#include "include/units.h"
using std::vector;
extern bool op_verbose;
extern bool op_batch;
OpticalPumping_Method::OpticalPumping_Method() {}
OpticalPumping_Method::OpticalPumping_Method(Eigenvector_Helper set_eigen,
                                             Laser_data set_laser_fe,
                                             Laser_data set_laser_ge)
  : eigen(set_eigen),
    numEStates(eigen.atom.numEStates), numFStates(eigen.atom.numFStates),
    numGStates(eigen.atom.numGStates), tau(eigen.atom.tau),
    laser_ge(set_laser_ge), laser_fe(set_laser_fe), Fe2_Vector(numEStates, 0),

    MFe2_Vector(numEStates, 0), Ff2_Vector(numFStates, 0),
    MFf2_Vector(numFStates, 0), Fg2_Vector(numGStates, 0),
    MFg2_Vector(numGStates, 0), nu_E(numEStates, 0.0), nu_F(numFStates, 0.0),
    nu_G(numGStates, 0.0),
    a_eg(numEStates, vector<vector<double> >(numGStates,
                                             vector<double>(3, 0.0))),
    a_ef(numEStates, vector<vector<double> >(numFStates,
                                             vector<double>(3, 0.0))),
    gFactor_E(numEStates, 0.0), cPlus_E(numEStates, 0.0),
    cPlus_F(numFStates, 0.0), cPlus_G(numGStates, 0.0),
    cMins_E(numEStates, 0.0), cMins_F(numFStates, 0.0),
    cMins_G(numGStates, 0.0) {
  dm_status = new DM_container(numEStates, numFStates, numGStates);
  dm_derivs = new DM_container(numEStates, numFStates, numGStates);
  // printf("\nOpticalPumping_Method::OpticalPumping_Method(...)\n");
  if (op_verbose) {
    printf("G = %d F = %d E = %d\n", numGStates, numFStates, numEStates);
  }
  setup_quantum_numbers(eigen.atom);
  setup_raising();
  setup_lowering();
  setup_gFactors(eigen.atom);
  setup_frequencies_excited(eigen.atom, eigen.field);
  setup_frequencies_ground(eigen.atom, eigen.field);
  setup_eg_coupling(eigen.atom);
  setup_ef_coupling(eigen.atom);
  if (op_verbose) print_couplings(stdout);
  setup_pop_uniform_ground();

  data.numGStates = eigen.atom.numGStates;
  data.numFStates = eigen.atom.numFStates;
  data.numEStates = eigen.atom.numEStates;
  data.atom = eigen.atom;
  data.laser_ge = set_laser_ge;
  data.laser_fe = set_laser_fe;
  data.a_eg = a_eg;
  data.a_ef = a_ef;
  data.nu_E = nu_E;
  data.nu_F = nu_F;
  data.nu_G = nu_G;
  // print_density_matrix(stdout);
  // print_data(stdout, 0.0);
  // for (int e = 0; e < numEStates; e++) printf("%15.10G\n",nu_E[e]);
}

OpticalPumping_Method::~OpticalPumping_Method() {}

void OpticalPumping_Method::update_population_euler(double dt) {
  reset_dPop();
  calculate_derivs(this -> dm_status);
  DM_container::mul(dm_derivs, dt);
  DM_container::add(dm_status, dm_derivs);
}

void OpticalPumping_Method::update_population_RK4(double dt) {

  calculate_derivs(this -> dm_status);
  DM_container *k1 = new DM_container(*dm_derivs);
  // printf("k1 = %g + %g i\n", GSL_REAL(k1->ff[0][0]), GSL_IMAG(k1->ff[0][0]));

  DM_container *arg = new DM_container(*k1);
  DM_container::mul(arg, dt/2.0);
  DM_container::add(arg, dm_status);
  calculate_derivs(arg);
  DM_container *k2 = new DM_container(*dm_derivs);

  arg = new DM_container(*k2);
  DM_container::mul(arg, dt/2.0);
  DM_container::add(arg, dm_status);
  calculate_derivs(arg);
  DM_container *k3 = new DM_container(*dm_derivs);

  arg = new DM_container(*k3);
  DM_container::mul(arg, dt);
  DM_container::add(arg, dm_status);
  calculate_derivs(arg);
  DM_container *k4 = new DM_container(*dm_derivs);

  DM_container *inc = new DM_container(dm_status->numEStates,
                                       dm_status->numFStates,
                                       dm_status->numGStates);
  DM_container::mul(k2, 2.0);
  DM_container::mul(k3, 2.0);
  
  DM_container::add(inc, k1);
  DM_container::add(inc, k2);
  DM_container::add(inc, k3);
  DM_container::add(inc, k4);
  DM_container::mul(inc, dt/6.0);
  DM_container::add(dm_status, inc);

  delete k1;
  delete k2;
  delete k3;
  delete k4;
  delete arg;
  delete inc;
}

void OpticalPumping_Method::setup_quantum_numbers(int I2, int Je2) {
  bool debug = false;
  int Fe2 = abs(I2 - Je2);
  int MFe2 = -Fe2;
  for (int e = 0; e < numEStates; e++) {
    Fe2_Vector[e] = Fe2;
    MFe2_Vector[e] = MFe2;
    if (debug) printf("e = %i \t Fe2 = %i \t Mfe2 = %i\n", e, Fe2_Vector[e],
           MFe2_Vector[e]);
    MFe2 += 2;
    if (MFe2 > Fe2) {
      Fe2 += 2;
      MFe2 = -Fe2;
    }
  }
  int Fg2 = (I2 - 1);
  int MFg2 = -Fg2;
  for (int g = 0; g < numGStates; g++) {
    Fg2_Vector[g] = Fg2;
    MFg2_Vector[g] = MFg2;
        if (debug) printf("g = %i \t Fg2 = %i \t Mfg2 = %i\n", g, Fg2_Vector[g],
         MFg2_Vector[g]);
    MFg2 += 2;
  }
  int Ff2 = (I2 + 1);
  int MFf2 = -Ff2;
  for (int f = 0; f < numFStates; f++) {
    Ff2_Vector[f] = Ff2;
    MFf2_Vector[f] = MFf2;
    if (debug) printf("f = %i \t Ff2 = %i \t Mff2 = %i\n", f, Ff2_Vector[f],
                  MFf2_Vector[f]);
    MFf2 += 2;
  }
}

void OpticalPumping_Method::setup_quantum_numbers(atom_data atom) {
  setup_quantum_numbers(atom.I2, atom.Je2);
}

void OpticalPumping_Method::setup_gFactors(atom_data atom) {
  Eigenvector_Helper alk;
  gFactor_G = alk.calc_gf(Fg2_Vector[0], 1, atom.I2, 0, 1, atom.g_I);
  gFactor_F = alk.calc_gf(Ff2_Vector[0], 1, atom.I2, 0, 1, atom.g_I);
  // printf("F-state g-Factor: %g\n", gFactor_G);
  // printf("G-state g-Factor: %g\n", gFactor_F);
  for (int e = 0; e < numEStates; e++) {
    gFactor_E[e] = alk.calc_gf(Fe2_Vector[e], atom.Je2, atom.I2, 2, 1,
                               atom.g_I);
    // printf("E-state[%d] g-Factor %g\n", e, gFactor_E[e]);
  }
}

void OpticalPumping_Method::setup_frequencies_excited(int I2, int Je2,
                                                      double excitation,
                                                      double hyperfine_const,
                                                      double B_z) {
  for (int e = 0; e < numEStates; e++) {
    nu_E[e] = set_frequency(excitation, I2, Je2, Fe2_Vector[e], MFe2_Vector[e],
                            hyperfine_const, B_z, gFactor_E[e]);
  }
}

void OpticalPumping_Method::setup_frequencies_excited(
                               atom_data atom,
                               magnetic_field_data field) {
  setup_frequencies_excited(atom.I2, atom.Je2, atom.nu_excited, atom.Aj_e,
                            field.B_z);
}

void OpticalPumping_Method::setup_frequencies_ground(int I2,
                                                     double hyperfine_const,
                                                     double B_z) {
  for (int g = 0; g < numGStates; g++) {
    nu_G[g] = set_frequency(0.0, I2, 1, Fg2_Vector[g], MFg2_Vector[g],
                            hyperfine_const, B_z, gFactor_G);
  }
  for (int f = 0; f < numFStates; f++) {
    nu_F[f] = set_frequency(0.0, I2, 1, Ff2_Vector[f], MFf2_Vector[f],
                            hyperfine_const, B_z, gFactor_F);
  }
}

void OpticalPumping_Method::setup_frequencies_ground(
                                atom_data atom,
                                magnetic_field_data field) {
  setup_frequencies_ground(atom.I2, atom.Aj_g, field.B_z);
}

double OpticalPumping_Method::set_frequency(double excitation, int I2, int J2,
                                            int F2, int Mf2,
                                            double hyperfine_const,
                                            double B_z, double g_f) {
  bool debug = false;
  // See DM equation 3.7
  double hyperfine = static_cast<double>((F2*(F2+2)) - (I2*(I2+2)) -
                                         (J2*(J2+2)));
  hyperfine *= hyperfine_const/8.0;
  if (debug) printf("Hyperfine: %15.10G MHz \t", hyperfine/_MHz);
  if (debug) printf("g_f = %6.4G, mu_B = %6.4G MHz/G, B_z = %6.4G G, ", g_f,
                    (_bohr_magneton/_planck_h)/(_MHz/_G), B_z/_G);
  double zeeman = static_cast<double>(Mf2);
  zeeman *= g_f * (_bohr_magneton/_planck_h) * B_z / 2.0;
  if (debug) printf("Zeeman: %15.10G MHz\t", zeeman/_MHz);
  if (debug) printf("Total: %15.10G MHz\n",
                    (excitation + hyperfine + zeeman)/_MHz);
  return (excitation + hyperfine + zeeman);
}

double OpticalPumping_Method::set_frequency(double excitation, int I2, int J2,
                                            int F2, int Mf2, int L2,
                                            double hyperfine_const,
                                            double g_I, double B_z) {
  Eigenvector_Helper alk;
  double g_f = alk.calc_gf(F2, J2, I2, L2, 1, g_I);
  return set_frequency(excitation, I2, J2, F2, Mf2, hyperfine_const, B_z, g_f);
}

void OpticalPumping_Method::setup_eg_coupling(int I2, int Je2) {
  for (int iq = 0; iq < 3; iq++) {
    int q2 = 2*(iq-1);  // iq holds the index in the array 1-3 and q holds the
    // actual polarization of the light:  sigma^- = -1; pi = 0; sigma^+ = +1
    for (int g = 0; g < numGStates; g++) {
      for (int e = 0; e < numEStates; e++) {
        a_eg[e][g][iq] = set_coupling(I2, 1, Fg2_Vector[g], MFg2_Vector[g], Je2,
                                     Fe2_Vector[e], MFe2_Vector[e], q2);
      }
    }
  }
}

void OpticalPumping_Method::setup_eg_coupling(atom_data atom) {
  setup_eg_coupling(atom.I2, atom.Je2);
}

void OpticalPumping_Method::setup_ef_coupling(int I2, int Je2) {
  for (int iq = 0; iq < 3; iq++) {
    int q2 = 2*(iq-1);             // See comment in setup_eg_couplings
    for (int f = 0; f < numFStates; f++) {
      for (int e = 0; e < numEStates; e++) {
        a_ef[e][f][iq] = set_coupling(I2, 1, Ff2_Vector[f], MFf2_Vector[f], Je2,
                                     Fe2_Vector[e], MFe2_Vector[e], q2);
      }
    }
  }
}

void OpticalPumping_Method::setup_ef_coupling(atom_data atom) {
  setup_ef_coupling(atom.I2, atom.Je2);
}

double OpticalPumping_Method::set_coupling(int I2, int Jg2, int Fg2, int Mg2,
                                           int Je2, int Fe2, int Me2, int q) {
  // The "q" parameter represents the true polarization of the light
  // -1->sigma^- 0->pi +1->sigma^+
  double expon = static_cast<double>(2+I2+Je2+Fe2+Fg2-Me2);
  expon /= 2.0;
  double coupling = pow(-1.0, expon);
  double temp = static_cast<double>((Fe2 + 1) * (Fg2 + 1) * (Je2 + 1));
  coupling *= sqrt(temp);
  coupling *= gsl_sf_coupling_3j(Fe2, 2, Fg2, -Me2, q, Mg2);
  coupling *= gsl_sf_coupling_6j(Fe2, 2, Fg2, Jg2, I2, Je2);
  return coupling;
}

void OpticalPumping_Method::setup_raising() {
  for (int e = 0; e < numEStates; e++) {
    cPlus_E[e] = set_raising(Fe2_Vector[e], MFe2_Vector[e]);
  }
  for (int f = 0; f < numFStates; f++) {
    cPlus_F[f] = set_raising(Ff2_Vector[f], MFf2_Vector[f]);
  }
  for (int g = 0; g < numGStates; g++) {
    cPlus_G[g] = set_raising(Fg2_Vector[g], MFg2_Vector[g]);
  }
}

void OpticalPumping_Method::setup_lowering() {
  for (int e = 0; e < numEStates; e++) {
    cMins_E[e] = set_lowering(Fe2_Vector[e], MFe2_Vector[e]);
  }
  for (int f = 0; f < numFStates; f++) {
    cMins_F[f] = set_lowering(Ff2_Vector[f], MFf2_Vector[f]);
  }
  for (int g = 0; g < numGStates; g++) {
    cMins_G[g] = set_lowering(Fg2_Vector[g], MFg2_Vector[g]);
  }
}

double OpticalPumping_Method::set_raising(int F2, int M2) {
  double F = (static_cast<double>(F2))/2.0;
  double M = (static_cast<double>(M2))/2.0;
  return sqrt((F-M)*(F+M+1));
}

double OpticalPumping_Method::set_lowering(int F2, int M2) {
  double F = (static_cast<double>(F2))/2.0;
  double M = (static_cast<double>(M2))/2.0;
  return sqrt((F+M)*(F-M+1));
}

void OpticalPumping_Method::setup_pop_uniform_ground() {
  double startPop = numFStates + numGStates;
  startPop = 1.0 / startPop;
  for (int f = 0; f < numFStates; f++) GSL_SET_REAL(&dm_status->ff[f][f],
                                                    startPop);
  for (int g = 0; g < numGStates; g++) GSL_SET_REAL(&dm_status->gg[g][g],
                                                    startPop);
}

void OpticalPumping_Method::print_couplings(FILE * des) {
  // First get minimum coupling to normalize everythign to integers
  // Then loop through each polarization, ground, excited state and print only
  // the couplings that are non-zero
  double min = 1.1;
  for (int q = 0; q < 3; q++) {
    for (int e = 0; e < numEStates; e++) {
      for (int f = 0; f < numFStates; f++) {
        if (fabs(a_ef[e][f][q]) < min && fabs(a_ef[e][f][q] > pow(10, -12)))
            min = fabs(a_ef[e][f][q]);
      }
      for (int g = 0; g < numGStates; g++) {
        if (fabs(a_eg[e][g][q]) < min && fabs(a_eg[e][g][q] > pow(10, -12)))
            min = fabs(a_eg[e][g][q]);
      }
    }
  }
  double norm = 1.0 / min;
  for (int q = 0; q < 3; q++) {
    fprintf(des, "**** q = %d ****\n", q-1);
    for (int g = 0; g < numGStates; g++) {
      for (int e = 0; e < numEStates; e++) {
        if (fabs(a_eg[e][g][q]) > pow(10, -6)) {
          fprintf(des, "|%d/2, %2d/2 --> |%d/2, %2d/2 >\t a = %8.6G\t",
                  Fg2_Vector[g], MFg2_Vector[g], Fe2_Vector[e], MFe2_Vector[e],
                  a_eg[e][g][q]);
          fprintf(des, "normalized = %8.6G\n", pow(norm*a_eg[e][g][q], 2.0));
        }
      }
    }
    for (int f = 0; f < numFStates; f++) {
      for (int e = 0; e < numEStates; e++) {
        if (fabs(a_ef[e][f][q]) > pow(10, -6)) {
          fprintf(des, "|%d/2, %2d/2 --> |%d/2, %2d/2 >\t a = %8.6G\t",
                  Ff2_Vector[f], MFf2_Vector[f], Fe2_Vector[e], MFe2_Vector[e],
                  a_ef[e][f][q]);
          fprintf(des, "normalized = %8.6G\n", pow(norm*a_ef[e][f][q], 2.0));
        }
      }
    }
  }
}

double OpticalPumping_Method::get_total_population() {
  double sum = 0.0;
  for (int g = 0; g < numGStates; g++) sum += GSL_REAL(dm_status->gg[g][g]);
  for (int f = 0; f < numFStates; f++) sum += GSL_REAL(dm_status->ff[f][f]);
  for (int e = 0; e < numEStates; e++) sum += GSL_REAL(dm_status->ee[e][e]);
  return sum;
}

double OpticalPumping_Method::get_polarization() {
  bool debug = false;
  // As in Pathria and Beale Ch5: <G> = Tr(rho*G) where G is any operator and
  // rho is the density matrix.  (They are matrix multiplied.)  For the nuclear
  // polarization, I use equation 3.25 from Dan's thesis.  Since I_z is
  // diagonal, I only need wory about the diagonal elements of the density
  // matrix.  In that case, I will forgo the matrix multiplication and just do
  // the sipmler method of sum of F-state{ sum over I,J eigenstates { I of that
  // eigenstate * mixing of F eigenstate into tha I,J eignestate * population of
  // that eignestate (diagonal density matrix element) } }
  int state;
  double state_pol = 0.0;
  for (state = 0; state < numGStates; state++) {
    if (debug) printf("For state %d, the decomp is:\n", state);
    for (int decomp = 0; decomp < eigen.atom.numBasisStates_ground; decomp++) {
      if (debug) {
        printf("\tI2 = %d --> %8.6G \n", eigen.nuclear_spin_ground[2*decomp],
               eigen.IzJz_decomp_ground[decomp][state]);
      }
      double tmp = (static_cast<double>(eigen.nuclear_spin_ground[2*decomp]));
      tmp /= 2.0;
      tmp *= pow(eigen.IzJz_decomp_ground[decomp][state], 2.0);
      tmp *= GSL_REAL(dm_status->gg[state][state]);
      state_pol += tmp;
    }
  }
  for (state = state; state < eigen.atom.numBasisStates_ground; state++) {
    if (debug) printf("For state %d, the decomp is:\n", state);
    for (int decomp = 0; decomp < eigen.atom.numBasisStates_ground; decomp++) {
      if (debug) {
        printf("\tI2 = %d --> %8.6G \n", eigen.nuclear_spin_ground[2*decomp],
               eigen.IzJz_decomp_ground[decomp][state]);
      }
      double tmp = (static_cast<double>(eigen.nuclear_spin_ground[2*decomp]));
      tmp /= 2.0;
      tmp *= pow(eigen.IzJz_decomp_ground[decomp][state], 2.0);
      tmp *= GSL_REAL(dm_status->ff[state-numGStates][state-numGStates]);
      state_pol += tmp;
    }
  }
  for (state = 0; state < eigen.atom.numBasisStates_excited; state++) {
    if (debug) printf("For e-state %d, the decomp is:\n", state);
    for (int decomp = 0; decomp < eigen.atom.numBasisStates_excited; decomp++) {
      if (debug) {
        printf("\tI2 = %d --> %8.6G \n", eigen.nuclear_spin_excited[2*decomp],
               eigen.IzJz_decomp_excited[decomp][state]);
      }
      double tmp = (static_cast<double>(eigen.nuclear_spin_excited[2*decomp]));
      tmp /= 2.0;
      tmp *= pow(eigen.IzJz_decomp_excited[decomp][state], 2.0);
      tmp *= GSL_REAL(dm_status->ee[state][state]);
      state_pol += tmp;
    }
  }
  state_pol /= ((static_cast<double>(eigen.atom.I2)) / 2.0);
  return state_pol;
}

double OpticalPumping_Method::get_alignment() {
  bool debug = false;
  // As in Pathria and Beale Ch5: <G> = Tr(rho*G) where G is any operator and
  // rho is the density matrix.  (They are matrix multiplied.)  For the nuclear
  // polarization, I use equation 3.26 from Dan's thesis.  Since I_z is
  // diagonal, I only need wory about the diagonal elements of the density
  // matrix.  In that case, I will forgo the matrix multiplication and just do
  // the sipmler method of sum of F-state{ sum over I,J eigenstates { I of that
  // eigenstate * mixing of F eigenstate into tha I,J eignestate * population of
  // that eignestate (diagonal density matrix element) } }
  int state;
  double state_align = 0.0;
  if (eigen.atom.I2 == 1) {      // My formula would divide by zero
    return 0.0;
  }
  for (state = 0; state < numGStates; state++) {
    if (debug) printf("For state %d, the decomp is:\n", state);
    for (int decomp = 0; decomp < eigen.atom.numBasisStates_ground; decomp++) {
      if (debug) {
        printf("\tI2 = %d --> %8.6G \n", eigen.nuclear_spin_ground[2*decomp],
               eigen.IzJz_decomp_ground[decomp][state]);
      }
      double tmp = (static_cast<double>(eigen.nuclear_spin_ground[2*decomp]));
      tmp /= 2.0;
      tmp = pow(tmp, 2.0);
      tmp *= pow(eigen.IzJz_decomp_ground[decomp][state], 2.0);
      tmp *= GSL_REAL(dm_status->gg[state][state]);
      state_align += tmp;
    }
  }
  for (state = state; state < eigen.atom.numBasisStates_ground; state++) {
    if (debug) printf("For state %d, the decomp is:\n", state);
    for (int decomp = 0; decomp < eigen.atom.numBasisStates_ground; decomp++) {
      if (debug) {
        printf("\tI2 = %d --> %8.6G \n", eigen.nuclear_spin_ground[2*decomp],
               eigen.IzJz_decomp_ground[decomp][state]);
      }
      double tmp = (static_cast<double>(eigen.nuclear_spin_ground[2*decomp]));
      tmp /= 2.0;
      tmp = pow(tmp, 2.0);
      tmp *= pow(eigen.IzJz_decomp_ground[decomp][state], 2.0);
      tmp *= GSL_REAL(dm_status->ff[state-numGStates][state-numGStates]);
      state_align += tmp;
    }
  }
  for (state = 0; state < eigen.atom.numBasisStates_excited; state++) {
    if (debug) printf("For e-state %d, the decomp is:\n", state);
    for (int decomp = 0; decomp < eigen.atom.numBasisStates_excited; decomp++) {
      if (debug) {
        printf("\tI2 = %d --> %8.6G \n", eigen.nuclear_spin_excited[2*decomp],
               eigen.IzJz_decomp_excited[decomp][state]);
      }
      double tmp = (static_cast<double>(eigen.nuclear_spin_excited[2*decomp]));
      tmp /= 2.0;
      tmp = pow(tmp, 2.0);
      tmp *= pow(eigen.IzJz_decomp_excited[decomp][state], 2.0);
      tmp *= GSL_REAL(dm_status->ee[state][state]);
      state_align += tmp;
    }
  }
  double nuc_spin = ((static_cast<double>(eigen.atom.I2)) / 2.0);
  state_align = nuc_spin*(nuc_spin + 1.0) - (3*state_align);
  state_align /= (nuc_spin*((2*nuc_spin) - 1.0));
  return state_align;
}

double OpticalPumping_Method::get_excited_state_total() {
  // Returns the sum of the excited state population.  This basically represents
  // the number of photions (unnormalized of course).
  double ePop = 0.0;
  for (int e = 0; e < numEStates; e++) {
    ePop += GSL_REAL(dm_status->ee[e][e]);
  }
  return ePop;
}

bool OpticalPumping_Method::is_hermitian() {
  // The density matrix is supposed to always remain Hermitian.   A hermitian
  // matrix has complex_conjugate(transpose(A)) == A.  What this means is the
  // real part must be symmetric and the imaginary part anti-symmetric.
  // Furthermore, the diagonal elements must be strictly real or the imaginary
  // part could not be anti-symmetric!
  bool hermit = true;
  double eps = pow(10, -8);     // Tolerance for hermitian-checking
  // First check that the diagonals are real
  for (int g = 0; g < numGStates; g++) {
    if (fabs(GSL_IMAG(dm_status->gg[g][g])) > eps) hermit = false;
  }
  for (int f = 0; f < numFStates; f++) {
    if (fabs(GSL_IMAG(dm_status->ff[f][f])) > eps) hermit = false;
  }
  for (int e = 0; e < numEStates; e++) {
    if (fabs(GSL_IMAG(dm_status->ee[e][e])) > eps) hermit = false;
  }

  // Now check that the off-diagonals have the appropriate symmetry
  for (int g = 0; g < numGStates; g++) {
    for (int gp = 0; gp < numGStates; gp++) {
      double imaginary_sum = fabs(GSL_IMAG(dm_status->gg[g][gp]) +
                                   GSL_IMAG(dm_status->gg[gp][g]));
      double real_diff = fabs(GSL_REAL(dm_status->gg[g][gp]) -
                               GSL_REAL(dm_status->gg[gp][g]));
      if ((imaginary_sum > eps) || (real_diff > eps)) hermit = false;
    }
  }
  for (int f = 0; f < numFStates; f++) {
    for (int fp = 0; fp < numFStates; fp++) {
      double imaginary_sum = fabs(GSL_IMAG(dm_status->ff[f][fp]) +
                                   GSL_IMAG(dm_status->ff[fp][f]));
      double real_diff = fabs(GSL_REAL(dm_status->ff[f][fp]) -
                               GSL_REAL(dm_status->ff[fp][f]));
      if ((imaginary_sum > eps) || (real_diff > eps)) hermit = false;
    }
  }
  for (int e = 0; e < numEStates; e++) {
    for (int ep = 0; ep < numGStates; ep++) {
      double imaginary_sum = fabs(GSL_IMAG(dm_status->ee[e][ep]) +
                                   GSL_IMAG(dm_status->ee[ep][e]));
      double real_diff = fabs(GSL_REAL(dm_status->ee[e][ep]) -
                               GSL_REAL(dm_status->ee[ep][e]));
      if ((imaginary_sum > eps) || (real_diff > eps)) hermit = false;
    }
  }
  return hermit;
}

void OpticalPumping_Method::print_data(FILE *des, double time) {
  double pop = get_total_population();
  double pol = get_polarization();
  double ali = get_alignment();
  double exc = get_excited_state_total();


  if (op_batch) {
    fprintf(des, "%8.6g\t", time/_us);
    fprintf(des, "%26.24G\t%26.24G\t%26.24G", exc, pol, ali);
  } else {
    fprintf(des, "%8.6G\t", time/_us);
    for (int g = 0; g < numGStates; g++) fprintf(des, "%8.6G\t",
                                                 GSL_REAL(dm_status->gg[g][g]));
    for (int f = 0; f < numFStates; f++) fprintf(des, "%8.6G\t",
                                                 GSL_REAL(dm_status->ff[f][f]));
    for (int e = 0; e < numEStates; e++) fprintf(des, "%8.6G\t",
                                                 GSL_REAL(dm_status->ee[e][e]));
    fprintf(des, "%8.6G\t%8.6G\t%8.6G\t%8.6G", pop, pol, ali, exc);
  }
    fprintf(des, "\n");
}

void OpticalPumping_Method::print_density_matrix(FILE *des) {
  // This function is useful when printing out the whole density matrix at one
  // time step to the screen.  It prints the real part as one matrix, tabs over
  // and prints the imaginary part as another matrix to the right of them first}
  fprintf(des, "rho_ee\n");
  for (int r = 0; r < numEStates; r++) {
    for (int c = 0; c < numEStates; c++) {
      fprintf(des, "%10.6G   ", GSL_REAL(dm_status->ee[r][c]));
    }
    fprintf(des, "\t\t");
    for (int c = 0; c < numEStates; c++) {
      fprintf(des, "%10.6G   ", GSL_IMAG(dm_status->ee[r][c]));
    }
    fprintf(des, "\n");
  }
  fprintf(des, "rho_gg\n");
  for (int r = 0; r < numGStates; r++) {
    for (int c = 0; c < numGStates; c++) {
      fprintf(des, "%10.6G   ", GSL_REAL(dm_status->gg[r][c]));
    }
    fprintf(des, "\t\t");
    for (int c = 0; c < numGStates; c++) {
      fprintf(des, "%10.6G   ", GSL_IMAG(dm_status->gg[r][c]));
    }
    fprintf(des, "\n");
  }
  fprintf(des, "rho_ff\n");
  for (int r = 0; r < numFStates; r++) {
    for (int c = 0; c < numFStates; c++) {
      fprintf(des, "%10.6G   ", GSL_REAL(dm_status->ff[r][c]));
    }
    fprintf(des, "\t\t");
    for (int c = 0; c < numFStates; c++) {
      fprintf(des, "%10.6G   ", GSL_IMAG(dm_status->ff[r][c]));
    }
    fprintf(des, "\n");
  }
  fprintf(des, "delta_rho_eg\n");
  for (int r = 0; r < numEStates; r++) {
    for (int c = 0; c < numGStates; c++) {
      fprintf(des, "%10.6G   ", GSL_REAL(dm_status->eg[r][c]));
    }
    fprintf(des, "\t\t");
    for (int c = 0; c < numGStates; c++) {
      fprintf(des, "%10.6G   ", GSL_IMAG(dm_status->eg[r][c]));
    }
    fprintf(des, "\n");
  }
  fprintf(des, "delta_rho_ef\n");
  for (int r = 0; r < numEStates; r++) {
    for (int c = 0; c < numFStates; c++) {
      fprintf(des, "%10.6G   ", GSL_REAL(dm_status->ef[r][c]));
    }
    fprintf(des, "\t\t");
    for (int c = 0; c < numFStates; c++) {
      fprintf(des, "%10.6G   ", GSL_IMAG(dm_status->ef[r][c]));
    }
    fprintf(des, "\n");
  }
  fprintf(des, "delta_fg\n");
  for (int r = 0; r < numFStates; r++) {
    for (int c = 0; c < numGStates; c++) {
      fprintf(des, "%10.6G   ", GSL_REAL(dm_status->fg[r][c]));
    }
    fprintf(des, "\t\t");
    for (int c = 0; c < numGStates; c++) {
      fprintf(des, "%10.6G   ", GSL_IMAG(dm_status->fg[r][c]));
    }
    fprintf(des, "\n");
  }
}

void OpticalPumping_Method::switch_off_laser(int las) {
  printf("Turning off laser %d\n", las);
  if (las == 1) {               // g --> e
    laser_ge.switch_off(tau);   // Need tau for saturation intensity
  } else if (las == 2) {        // f --> e
    laser_fe.switch_off(tau);   // Need tau for saturation intensity
  } else {
    printf("***Warning***\nAttmpeting to switch off laser %d.", las);
    printf(", but there is no such laser.  No action taken\n***Warning***\n");
  }
}

void OpticalPumping_Method::change_magnetic_field(double newField) {
  // printf("OP::change_magnetic_field\n");
  eigen.field.B_z = newField;
  setup_frequencies_excited(eigen.atom, eigen.field);
  setup_frequencies_ground(eigen.atom, eigen.field);
  eigen.IzJz_decomp_ground = eigen.diagH(0);
  eigen.IzJz_decomp_excited = eigen.diagH(1);
}

void OpticalPumping_Method::reset_dPop() {
  for (int e = 0; e < numEStates; e++) {
    for (int ep = 0; ep < numEStates; ep++) {
      dm_derivs->ee[e][ep] = gsl_complex_rect(0.0, 0.0);
    }
    for (int f = 0; f < numFStates; f++) {
      dm_derivs->ef[e][f] = gsl_complex_rect(0.0, 0.0);
    }
    for (int g = 0; g < numGStates; g++) {
      dm_derivs->eg[e][g] = gsl_complex_rect(0.0, 0.0);
    }
  }
  for (int f = 0; f < numFStates; f++) {
    for (int fp = 0; fp < numFStates; fp++) {
      dm_derivs->ff[f][fp] = gsl_complex_rect(0.0, 0.0);
    }
    for (int g = 0; g < numGStates; g++) {
      dm_derivs->fg[f][g] = gsl_complex_rect(0.0, 0.0);
    }
  }
  for (int g = 0; g < numGStates; g++) {
    for (int gp = 0; gp < numGStates; gp++) {
      dm_derivs->gg[g][gp] = gsl_complex_rect(0.0, 0.0);
    }
  }
}

void OpticalPumping_Method::apply_dPop(double dt) {
  for (int e = 0; e < numEStates; e++) {
    for (int ep = 0; ep < numEStates; ep++) {
      dm_derivs->ee[e][ep] = gsl_complex_mul_real(dm_derivs->ee[e][ep], dt);
      dm_status->ee[e][ep] = gsl_complex_add(dm_status->ee[e][ep],
                                             dm_derivs->ee[e][ep]);
      if ((e == ep) && GSL_REAL(dm_status->ee[e][ep]) < 0.0) {
        dm_status->ee[e][ep] = gsl_complex_rect(0.0, 0.0);
      }
    }
    for (int f = 0; f < numFStates; f++) {
      dm_derivs->ef[e][f] = gsl_complex_mul_real(dm_derivs->ef[e][f], dt);
      dm_status->ef[e][f] = gsl_complex_add(dm_status->ef[e][f],
                                            dm_derivs->ef[e][f]);
    }
    for (int g = 0; g < numGStates; g++) {
      dm_derivs->eg[e][g] = gsl_complex_mul_real(dm_derivs->eg[e][g], dt);
      dm_status->eg[e][g] = gsl_complex_add(dm_status->eg[e][g],
                                            dm_derivs->eg[e][g]);
    }
  }
  for (int f = 0; f < numFStates; f++) {
    for (int fp = 0; fp < numFStates; fp++) {
      dm_derivs->ff[f][fp] = gsl_complex_mul_real(dm_derivs->ff[f][fp], dt);
      dm_status->ff[f][fp] = gsl_complex_add(dm_status->ff[f][fp],
                                             dm_derivs->ff[f][fp]);
      if ((f == fp) && GSL_REAL(dm_status->ff[f][fp]) < 0.0) {
        dm_status->ff[f][fp] = gsl_complex_rect(0.0, 0.0);
      }
    }
    for (int g = 0; g < numGStates; g++) {
      dm_derivs->fg[f][g] = gsl_complex_mul_real(dm_derivs->fg[f][g], dt);
      dm_status->fg[f][g] = gsl_complex_add(dm_status->fg[f][g],
                                            dm_derivs->fg[f][g]);
    }
  }
  for (int g = 0; g < numGStates; g++) {
    for (int gp = 0; gp < numGStates; gp++) {
      dm_derivs->gg[g][gp] = gsl_complex_mul_real(dm_derivs->gg[g][gp], dt);
      dm_status->gg[g][gp] = gsl_complex_add(dm_status->gg[g][gp],
                                             dm_derivs->gg[g][gp]);
      if ((g == gp) && GSL_REAL(dm_status->gg[g][gp]) < 0.0) {
        dm_status->gg[g][gp] = gsl_complex_rect(0.0, 0.0);
      }
    }
  }
}
