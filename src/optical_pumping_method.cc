// Authors: Benjamin Fenker 2013
// Copyright 2012 Benjamin Fenker

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <gsl/gsl_sf_coupling.h>
#include "include/optical_pumping_method.h"
#include "include/units.h"
using std::string;
using std::vector;
extern bool op_verbose;
extern bool op_batch;
OpticalPumping_Method::OpticalPumping_Method() {}
OpticalPumping_Method::OpticalPumping_Method(Eigenvector_Helper set_eigen,
                                             Laser_data set_laser_fe,
                                             Laser_data set_laser_ge,
                                             double tilt)
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
  // laser_fe.nu = (nu_E[5] - nu_F[2]) + laser_fe.detune;
  // laser_ge.nu = (nu_E[5] - nu_G[1]) + laser_ge.detune;
  // printf("Excited 22 energy: %8.6f MHz\n", nu_E[5]);
  // printf("Ground  20 energy: %8.6f MHz\n", nu_F[2]);
  // printf("Ground  10 energy: %8.6f MHz\n", nu_G[1]);
  // printf("Forced ge frequency to: %12.10g MHz\n", laser_ge.nu/_MHz);
  // printf("Forced fe frequency to: %12.10g MHz\n", laser_fe.nu/_MHz);
  setup_eg_coupling(eigen.atom);
  setup_ef_coupling(eigen.atom);
  if (op_verbose) print_couplings(stdout);
  // setup_pop_uniform_ground();
  //  setup_pop_withTilt(tilt);
  tilt = tilt;
  setup_pop_arbitrary("initial.pop");

  data.numGStates = eigen.atom.numGStates;
  data.numFStates = eigen.atom.numFStates;
  data.numEStates = eigen.atom.numEStates;
  data.atom = eigen.atom;
  data.laser_ge = set_laser_ge;
  data.laser_fe = set_laser_fe;
  // data.laser_ge = laser_ge;
  // data.laser_fe = laser_fe;
  data.a_eg = a_eg;
  data.a_ef = a_ef;
  data.nu_E = nu_E;
  data.nu_F = nu_F;
  data.nu_G = nu_G;
  // print_density_matrix(stdout);
  // print_data(stdout, 0.0);
  // for (int e = 0; e < numEStates; e++) printf("%15.10G\n",nu_E[e]);
}

OpticalPumping_Method::~OpticalPumping_Method() {
  // For some reason this crashes things on the second go-round
  // I'm not sure why.  But this is SUSPICIOUS!
  // if (dm_status) delete dm_status;
  // if (dm_derivs) delete dm_derivs;
}

void OpticalPumping_Method::update_population_euler(double dt) {
  double dt_orig = dt;
  reset_dPop();
  calculate_derivs(this -> dm_status);
  // printf("Derivative for stretched state is\n%g ",
  //        GSL_REAL(dm_derivs->ff[4][4]));  
  DM_container::mul(dm_derivs, dt);
  // while (!DM_container::okay_to_add(dm_status, dm_daerivs)) {
  //   printf("Problems at dt = %g ns, ", dt_orig/_ns);
  //   dt = dt*0.5;
  //   if (dt < 0.001 *_ns) exit(0);
  //   printf("trying dt = %g ns\n", dt/_ns);
  //   DM_container::mul(dm_derivs, 0.5);
  // }
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

  delete arg;
  arg = new DM_container(*k2);
  DM_container::mul(arg, dt/2.0);
  DM_container::add(arg, dm_status);
  calculate_derivs(arg);
  DM_container *k3 = new DM_container(*dm_derivs);

  delete arg;
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
  bool j = inc -> equalsZero();
  DM_container::mul(inc, dt/6.0);

  if (j) {
    // printf("Zero!\n");
  }
  //  if (inc -> equalsZero()) printf("Zero!\n");
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
  // printf("F-state g-Factor: %g\n", gFactor_F);
  // printf("G-state g-Factor: %g\n", gFactor_G);
  for (int e = 0; e < numEStates; e++) {
    gFactor_E[e] = alk.calc_gf(Fe2_Vector[e], atom.Je2, atom.I2, 2, 1,
                               atom.g_I);
    //       printf("E-state[%d] g-Factor %g\n", e, gFactor_E[e]);
  }
}

void OpticalPumping_Method::setup_frequencies_excited(int I2, int Je2,
                                                      double excitation,
                                                      double hyperfine_const,
                                                      double B_z) {
  //  printf("Excited state frequencies...\n");
  for (int e = 0; e < numEStates; e++) {
    nu_E[e] = set_frequency(excitation, I2, Je2, Fe2_Vector[e], MFe2_Vector[e],
                            hyperfine_const, B_z, gFactor_E[e]);
    // printf("|%g, %g > = %g _MHz\n", static_cast<double>(Fe2_Vector[e])/2.0,
    //        static_cast<double>(MFe2_Vector[e])/2.0, (nu_E[e] - excitation)/_MHz);

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
  //  printf("Ground state energies...\n");
  for (int g = 0; g < numGStates; g++) {
    nu_G[g] = set_frequency(0.0, I2, 1, Fg2_Vector[g], MFg2_Vector[g],
                            hyperfine_const, B_z, gFactor_G);
    // printf("|%g, %g > = %g _MHz\n", static_cast<double>(Fg2_Vector[g])/2.0,
    //        static_cast<double>(MFg2_Vector[g])/2.0, nu_G[g]/_MHz);
  }
  for (int f = 0; f < numFStates; f++) {
    nu_F[f] = set_frequency(0.0, I2, 1, Ff2_Vector[f], MFf2_Vector[f],
                            hyperfine_const, B_z, gFactor_F);
    // printf("|%g, %g > = %g _MHz\n", static_cast<double>(Ff2_Vector[f])/2.0,
    //        static_cast<double>(MFf2_Vector[f])/2.0, nu_F[f]/_MHz);

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
  bool debug = false;;
  // See DM equation 3.7
  double hyperfine = static_cast<double>((F2*(F2+2)) - (I2*(I2+2)) -
                                         (J2*(J2+2)));
  hyperfine *= hyperfine_const/8.0;
  if (debug) printf("I2 = %d\t J2 = %d\t F2 = %d\t Mf2 = %d\n",
                    I2, J2, F2, Mf2);
  if (debug) printf("Hyperfine: %6.4f MHz \t", hyperfine/_MHz);
  if (debug) printf("g_f = %6.4G, mu_B = %6.4G MHz/G, B_z = %6.4G G, ", g_f,
                    (_bohr_magneton/_planck_h)/(_MHz/_G), B_z/_G);
  double zeeman = static_cast<double>(Mf2);
  zeeman *= -1.0*g_f * (_bohr_magneton/_planck_h) * B_z / 2.0;
  if (debug) printf("Zeeman: %6.4f MHz\t", zeeman/_MHz);
  if (debug) printf("Total: %6.4f MHz\n",
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

void OpticalPumping_Method::setup_pop_arbitrary(string fname) {
  std::ifstream ifs;
  ifs.open(fname.c_str(), std::ifstream::in);
  if (!ifs.good()) {
    printf("Could not open population file %s\n", fname.c_str());
    setup_pop_uniform_ground();
    return;
  }
  vector<double> pop(numFStates+numGStates, 0.0);
  double sum = 0.0;
  for (int i = 0; i < numFStates + numGStates; i++) {
    if (!ifs.good()) {
      printf("Population file %s not enough inputs\n", fname.c_str());
      setup_pop_uniform_ground();
      return;
    }
    double tmp;
    ifs >> tmp;
    if (tmp >= 0) {
      pop[i] = tmp;
    } else {
      printf("WARNING: Initial population < 0 in file! ");
      printf("(Set this population to zero and continue)\n");
      pop[i] = 0;
    }
    sum += pop[i];
  }

  if (sum <= 0) {
    printf("ERROR... Total initial population <= 0... set to uniform");
    setup_pop_uniform_ground();
    return;
  }

  for (int i = 0; i < numFStates + numGStates; i++) pop[i] /= sum;

  //  printf("Initial populations...\n");
  for (int i = 0; i < numGStates; i++) {
    GSL_SET_REAL(&dm_status->gg[i][i], pop[i]);
    // printf("| F = %d/2, Mf = %d/2 > = %g\t", Fg2_Vector[i], MFg2_Vector[i],
    //        pop[i]);
  }
  printf("\n");
  for (int i = 0; i < numFStates; i++) {
    GSL_SET_REAL(&dm_status->ff[i][i], pop[i+numGStates]);
    // printf("| F = %d/2, Mf = %d/2 > = %g\t", Ff2_Vector[i], MFf2_Vector[i],
    //        pop[i+numGStates]);
  }
  // printf("\n");
}

void OpticalPumping_Method::setup_pop_withTilt(double tilt) {
  double startPop = numFStates + numGStates;
  startPop = 1.0 / startPop;
  // printf("First comparison: %g\n", static_cast<double>(numFStates)/2.0);
  for (int f = 0; f < floor(static_cast<double>(numFStates)/2.0); f++) {
    GSL_SET_REAL(&dm_status->ff[f][f], startPop*(1.0-tilt));
    // printf("Filling index %d with %g\n", f, (1.0-tilt));
  }
  int highStarting = numFStates/2;
  if (numFStates % 2 != 0) {
    GSL_SET_REAL(&dm_status->ff[(numFStates/2)][(numFStates/2)],
                 startPop);
    // printf("Filling index %d with %g\n", (numFStates/2), 1.0);
    highStarting++;
  }

  // printf("High starting: %d\n", highStarting);
  for (int f = highStarting; f < numFStates; f++) {
    GSL_SET_REAL(&dm_status->ff[f][f], startPop*(1.0+tilt));
    // printf("Filling index %d with %g\n", f, (1.0+asym));
  }

  for (int g = 0; g < floor(static_cast<double>(numGStates)/2.0); g++) {
    GSL_SET_REAL(&dm_status->gg[g][g], startPop*(1.0-tilt));
    // printf("Filling index %d with %g\n", g, (1.0-tilt));
  }
  highStarting = numGStates/2;
  if (numGStates % 2 != 0) {
    GSL_SET_REAL(&dm_status->gg[(numGStates/2)][(numGStates/2)],
                 startPop);
    // printf("Filling index %d with %g\n", (numGStates/2), 1.0);
    highStarting++;
  }

  // printf("High starting: %d\n", highStarting);
  for (int g = highStarting; g < numGStates; g++) {
    GSL_SET_REAL(&dm_status->gg[g][g], startPop*(1.0+tilt));
    // printf("Filling index %d with %g\n", g, (1.0+asym));
  }
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
    if (debug) printf("For g-state %d, the decomp is:\n", state);
    for (int decomp = 0; decomp < eigen.atom.numBasisStates_ground; decomp++) {
      if (debug) {
        printf("\tI2 = %d --> %8.6G \t", eigen.nuclear_spin_ground[2*decomp],
               eigen.IzJz_decomp_ground[decomp][state]);
      }
      double tmp = (static_cast<double>(eigen.nuclear_spin_ground[2*decomp]));
      tmp /= 2.0;
      tmp *= pow(eigen.IzJz_decomp_ground[decomp][state], 2.0);
      tmp *= GSL_REAL(dm_status->gg[state][state]);
      if (debug) printf("And population: %g, contributes %g\n",
                        GSL_REAL(dm_status->gg[state][state]), tmp);
      state_pol += tmp;
    }
  }
  for (state = state; state < eigen.atom.numBasisStates_ground; state++) {
    if (debug) printf("For f-state %d, the decomp is:\n", state);
    for (int decomp = 0; decomp < eigen.atom.numBasisStates_ground; decomp++) {
      if (debug) {
        printf("\tI2 = %d --> %8.6G \t", eigen.nuclear_spin_ground[2*decomp],
               eigen.IzJz_decomp_ground[decomp][state]);
      }
      double tmp = (static_cast<double>(eigen.nuclear_spin_ground[2*decomp]));
      tmp /= 2.0;
      tmp *= pow(eigen.IzJz_decomp_ground[decomp][state], 2.0);
      tmp *= GSL_REAL(dm_status->ff[state-numGStates][state-numGStates]);
      if (debug) {
        printf("And population: %g, contributes %g\n",
               GSL_REAL(dm_status->ff[state-numGStates][state-numGStates]),
               tmp);
      }

      state_pol += tmp;
    }
  }
  for (state = 0; state < eigen.atom.numBasisStates_excited; state++) {
    if (debug) printf("For e-state %d, the decomp is:\n", state);
    for (int decomp = 0; decomp < eigen.atom.numBasisStates_excited; decomp++) {
      if (debug) {
        // printf("\tI2 = %d --> %8.6G \t", eigen.nuclear_spin_excited[2*decomp],
        //        eigen.IzJz_decomp_excited[decomp][state]);
      }
      double tmp = (static_cast<double>(eigen.nuclear_spin_excited[2*decomp]));
      tmp /= 2.0;
      tmp *= pow(eigen.IzJz_decomp_excited[decomp][state], 2.0);
      tmp *= GSL_REAL(dm_status->ee[state][state]);
      // if (debug) printf("And population: %g, contributes %g\n",
      //                   GSL_REAL(dm_status->ee[state][state]), tmp);
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
  double eps = pow(10, -6);     // Tolerance for hermitian-checking
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

  // Now check that the populations are between 0 and 1
  for (int e = 0; e < numEStates; e++) {
    if (GSL_REAL(dm_status->ee[e][e]) < (0.0 - eps)
        || GSL_REAL(dm_status->ee[e][e]) > (1.0 + eps)) {
      hermit = false;
    }
  }
  for (int f = 0; f < numFStates; f++) {
    if (GSL_REAL(dm_status->ff[f][f]) < (0.0 - eps)
        || GSL_REAL(dm_status->ff[f][f]) > (1.0 + eps)) {
      hermit = false;
    }
  }
  for (int g = 0; g < numGStates; g++) {
    if (GSL_REAL(dm_status->gg[g][g]) < (0.0 - eps)
        || GSL_REAL(dm_status->gg[g][g]) > (1.0 + eps)) {
      hermit = false;
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
  printf("OP::change_magnetic_field\n");
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

void OpticalPumping_Method::apply_transverse_field_ee(DM_container *in) {
#pragma omp parallel for
  for (int e = 0; e < numEStates; e++) {
    for (int ep = 0; ep < numEStates; ep++) {
      // if ((e == ep) ||
      //     (es_hyperfine && MFe2_Vector[e] == MFe2_Vector[ep]) ||
      //     (es_Zeeman && Fe2_Vector[e] == Fe2_Vector[ep]) ||
      //     (es_hyperfine && es_Zeeman)) {
        gsl_complex Bx = gsl_complex_rect(0.0, 0.0);
        gsl_complex left = gsl_complex_rect(0.0, 0.0);
        gsl_complex rigt = gsl_complex_rect(0.0, 0.0);
        // printf("For e = %d, comparing %d to %d\n", e, MFe2_Vector[e],
        //        MFe2_Vector[e+1]);
        // Is there a state with Mf` = Mf+1 and F` = F?
        if (e+1 < numEStates) {         // Keep things in bounds
          if (MFe2_Vector[e]+2 == MFe2_Vector[e+1]) {  // Sanity check
            gsl_complex tmp = in->ee[e+1][ep];
            // printf("Bx = %g + %g i + ", GSL_REAL(tmp), GSL_IMAG(tmp));
            tmp = gsl_complex_mul_real(tmp, cPlus_E[e]);
            left = gsl_complex_add(left, tmp);
          }
        }
        // Is there a state with Mf` = Mf-1 and F` = F?
        if (e > 0) {                    // Keep things in bounds
          if (MFe2_Vector[e]-2 == MFe2_Vector[e-1]) {  // Sanity check
            gsl_complex tmp = in->ee[e-1][ep];
            // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
            tmp = gsl_complex_mul_real(tmp, cMins_E[e]);
            left = gsl_complex_add(left, tmp);
          }
        }
        left = gsl_complex_mul_real(left, gFactor_E[e]);
        // Is there a state with Mf` = Mf+1 and F` = F?
        if (ep+1 < numEStates) {         // Keep things in bounds
          if (MFe2_Vector[ep]+2 == MFe2_Vector[ep+1]) {  // Sanity check
            gsl_complex tmp = in->ee[e][ep+1];
            // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
            tmp = gsl_complex_mul_real(tmp, cPlus_E[ep]);
            rigt = gsl_complex_sub(rigt, tmp);
          }
        }
        // Is there a state with Mf` = Mf-1 and F` = F?
        if (ep > 0) {         // Keep things in bounds
          if (MFe2_Vector[ep]-2 == MFe2_Vector[ep-1]) {
            gsl_complex tmp = in->ee[e][ep-1];
            // printf("%g + %g i\n", GSL_REAL(tmp), GSL_IMAG(tmp));
            tmp = gsl_complex_mul_real(tmp, cMins_E[ep]);
            rigt = gsl_complex_sub(rigt, tmp);
          }
        }
        rigt = gsl_complex_mul_real(rigt, gFactor_E[ep]);
        Bx = gsl_complex_add(left, rigt);
        Bx = gsl_complex_mul_imag(Bx,
                                  -_bohr_magneton*eigen.field.B_x
                                  /2.0/_planck_hbar);
        // printf("Field is: %g\t", eigen.field.B_x/_G);
        // printf("Bx contribution is %g + %g i\n", GSL_REAL(Bx), GSL_IMAG(Bx));
        dm_derivs->ee[e][ep] = gsl_complex_add(dm_derivs->ee[e][ep], Bx);
      // }   // End if coherences
    }      // End ep
  }        // End e
}

void OpticalPumping_Method::apply_transverse_field_gg(DM_container *in) {
#pragma omp parallel for
  for (int g = 0; g < numGStates; g++) {
    for (int gp = 0; gp < numGStates; gp++) {
      if ((g == gp) || (gs_Zeeman)) {
        // Get ready for transverse magnetic field!
        gsl_complex Bx = gsl_complex_rect(0.0, 0.0);
        gsl_complex left = gsl_complex_rect(0.0, 0.0);
        gsl_complex rigt = gsl_complex_rect(0.0, 0.0);
        // printf("For e = %d, comparing %d to %d\n", e, MFe2_Vector[e],
        //        MFe2_Vector[e+1]);
        // Is there a state with Mf` = Mf+1 and F` = F?
        if (g+1 < numGStates) {         // Keep things in bounds
          if (MFg2_Vector[g]+2 == MFg2_Vector[g+1]) {  // Sanity check
            gsl_complex tmp = in->gg[g+1][gp];
            // printf("Bx = %g + %g i + ", GSL_REAL(tmp), GSL_IMAG(tmp));
            tmp = gsl_complex_mul_real(tmp, cPlus_G[g]);
            left = gsl_complex_add(left, tmp);
          }
        }
        // Is there a state with Mf` = Mf-1 and F` = F?
        if (g > 0) {         // Keep things in bounds
          if (MFg2_Vector[g]-2 == MFg2_Vector[g-1]) {  // Sanity check
            gsl_complex tmp = in->gg[g-1][gp];
            // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
            tmp = gsl_complex_mul_real(tmp, cMins_G[g]);
            left = gsl_complex_add(left, tmp);
          }
        }
        left = gsl_complex_mul_real(left, gFactor_G);
        // Is there a state with Mf` = Mf+1 and F` = F?
        if (gp + 1 < numGStates) {      // Keep things in bounds
          if (MFg2_Vector[gp]+2 == MFg2_Vector[gp+1]) {  // Sanity check
            gsl_complex tmp = in->gg[g][gp+1];
            // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
            tmp = gsl_complex_mul_real(tmp, cPlus_G[gp]);
            rigt = gsl_complex_sub(rigt, tmp);
          }
        }
        // Is there a state with Mf` = Mf-1 and F` = F?
        if (gp > 0) {                   // Keep things in bounds
          if (MFg2_Vector[gp]-2 == MFg2_Vector[gp-1]) {  // Sanity check
            gsl_complex tmp = in->gg[g][gp-1];
            // printf("%g + %g i\n", GSL_REAL(tmp), GSL_IMAG(tmp));
            tmp = gsl_complex_mul_real(tmp, cMins_G[gp]);
            rigt = gsl_complex_sub(rigt, tmp);
          }
        }
        rigt = gsl_complex_mul_real(rigt, gFactor_G);
        Bx = gsl_complex_add(left, rigt);
        Bx = gsl_complex_mul_imag(Bx,
                                  -_bohr_magneton*eigen.field.B_x
                                  /2.0/_planck_hbar);
        // printf("Bx contribution is %g + %g i\n", GSL_REAL(Bx), GSL_IMAG(Bx));
        dm_derivs->gg[g][gp] = gsl_complex_add(dm_derivs->gg[g][gp], Bx);
      }  // End if coherences
    }    // End gp
  }      // End g
}        // End B_x gg

/*
void OpticalPumping_Method::apply_transverse_field_ff(DM_container *in) {

  for (int i = 0; i < numFStates; i++) {
    for (int j = 0; j < numFStates; j++) {
      gsl_complex bx_contribution = gsl_complex_rect(0.0, 0.0);
      double c, mf;
      gsl_complex rho;

      mf = (static_cast<double>(MFf2_Vector[i]))/2.0;
      c = sqrt( (2-mf) * (2+mf+1) );

      //      printf ("For i = %d, j = %d \t c = %g\n", i, j, c);
      if (fabs(c) > 0.0) {
        rho = in -> ff[i+1][j];
        bx_contribution = gsl_complex_add(bx_contribution,
                                          gsl_complex_mul_real(rho, c));
      }

      c = sqrt( (2+mf) * (2-mf+1) );
      if (fabs(c) > 0.0) {
        rho = in -> ff[i-1][j];
        bx_contribution = gsl_complex_add(bx_contribution,
                                          gsl_complex_mul_real(rho, c));
      }

      mf = (static_cast<double>(MFf2_Vector[j]))/2.0;
      c = sqrt( (2-mf) * (2+mf+1) );
      if (fabs(c) > 0.0) {
        rho = in -> ff[i][j+1];
        bx_contribution = gsl_complex_sub(bx_contribution,
                                          gsl_complex_mul_real(rho, c));
      }

      c = sqrt( (2+mf) * (2-mf+1) );
      if (fabs(c) > 0.0) {
        rho = in -> ff[i][j-1];
        bx_contribution = gsl_complex_sub(bx_contribution,
                                          gsl_complex_mul_real(rho, c));
      }

      double mul_factor =
          -1.0*_bohr_magneton*eigen.field.B_x*gFactor_F / (2.0*_planck_hbar);
      bx_contribution = gsl_complex_mul_imag(bx_contribution, mul_factor);

      dm_derivs->ff[i][j] =
          gsl_complex_add(dm_derivs->ff[i][j], bx_contribution);
    }
  }
}
*/

/*
void OpticalPumping_Method::apply_transverse_field_ff(DM_container *in) {
  for (int f = 0; f < numFStates; f++) {
    for (int fp = 0; fp < numFStates; fp++) {
        // Get ready for transverse magnetic field!
        gsl_complex Bx = gsl_complex_rect(0.0, 0.0);
        // gsl_complex left = gsl_complex_rect(0.0, 0.0);
        // gsl_complex rigt = gsl_complex_rect(0.0, 0.0);
        // Is there a state with Mf` = Mf+1 and F` = F?
        double c, mf;

        mf = (static_cast<double>(MFf2_Vector[f]))/2.0;
        c = sqrt( (2-mf) * (2+mf+1) );
        if (fabs(c) > 0) {
            gsl_complex tmp = in->ff[f+1][fp];

            // Bx = gsl_complex_add(Bx, gsl_complex_mul_real(tmp, c));
            tmp = gsl_complex_mul_real(tmp, c);
             Bx = gsl_complex_add(Bx, tmp);
        }


        c = sqrt( (2+mf) * (2-mf+1) );        
        // Is there a state with Mf` = Mf-1 and F` = F?
        if (fabs(c) > 0) {                    // Keep things in bounds
            gsl_complex tmp = in->ff[f-1][fp];

            // Bx = gsl_complex_add(Bx, gsl_complex_mul_real(tmp, c));

            tmp = gsl_complex_mul_real(tmp, c);
            Bx = gsl_complex_add(Bx, tmp);
        }

        // Is there a state with Mf` = Mf+1 and F` = F?
        mf = (static_cast<double>(MFf2_Vector[fp]))/2.0;
        c = sqrt( (2-mf) * (2+mf+1) );

        if (fabs(c) > 0) {        // Keep things in bounds
            gsl_complex tmp = in->ff[f][fp+1];

            //Bx = gsl_complex_sub(Bx, gsl_complex_mul_real(tmp, c));
            tmp = gsl_complex_mul_real(tmp, c);
            Bx = gsl_complex_sub(Bx, tmp);
        }
        // Is there a state with Mf` = Mf-1 and F` = F?
        c = sqrt( (2+mf) * (2-mf+1) );        
        if (fabs(c) > 0) {                   // Keep things in bounds
            gsl_complex tmp = in->ff[f][fp-1];
            //     Bx = gsl_complex_sub(Bx, gsl_complex_mul_real(tmp, c));
            tmp = gsl_complex_mul_real(tmp, c);
            Bx = gsl_complex_sub(Bx, tmp);
        }

      double mul_factor =
          -1.0*_bohr_magneton*eigen.field.B_x*gFactor_F / (2.0*_planck_hbar);
      Bx = gsl_complex_mul_imag(Bx, mul_factor);

      dm_derivs->ff[f][fp] = gsl_complex_add(dm_derivs->ff[f][fp], Bx);
    }    // End fp
  }      // End f
}        // End function
*/


void OpticalPumping_Method::apply_transverse_field_ff(DM_container *in) {
#pragma omp parallel for
  for (int f = 0; f < numFStates; f++) {
    for (int fp = 0; fp < numFStates; fp++) {
        // Get ready for transverse magnetic field!
        gsl_complex Bx = gsl_complex_rect(0.0, 0.0);
        gsl_complex left = gsl_complex_rect(0.0, 0.0);
        gsl_complex rigt = gsl_complex_rect(0.0, 0.0);
        // Is there a state with Mf` = Mf+1 and F` = F?
        if (f+1 < numFStates) {         // Keep things in bounds
            gsl_complex tmp = in->ff[f+1][fp];
            tmp = gsl_complex_mul_real(tmp, cPlus_F[f]);
            left = gsl_complex_add(left, tmp);
        }
        // Is there a state with Mf` = Mf-1 and F` = F?
        if (f > 0) {                    // Keep things in bounds
            gsl_complex tmp = in->ff[f-1][fp];
            tmp = gsl_complex_mul_real(tmp, cMins_F[f]);
            left = gsl_complex_add(left, tmp);
        }
        left = gsl_complex_mul_real(left, gFactor_F);
        // Is there a state with Mf` = Mf+1 and F` = F?
        if (fp+1 < numFStates) {        // Keep things in bounds
            gsl_complex tmp = in->ff[f][fp+1];
            tmp = gsl_complex_mul_real(tmp, cPlus_F[fp]);
            rigt = gsl_complex_add(rigt, tmp);
        }
        // Is there a state with Mf` = Mf-1 and F` = F?
        if (fp > 0) {                   // Keep things in bounds
            gsl_complex tmp = in->ff[f][fp-1];
            tmp = gsl_complex_mul_real(tmp, cMins_F[fp]);
            rigt = gsl_complex_add(rigt, tmp);
        }
        rigt = gsl_complex_mul_real(rigt, gFactor_F);
        Bx = gsl_complex_sub(left, rigt);
        Bx = gsl_complex_mul_imag(Bx,
                                  -_bohr_magneton*eigen.field.B_x
                                  /2.0/_planck_hbar);

        dm_derivs->ff[f][fp] = gsl_complex_add(dm_derivs->ff[f][fp], Bx);
    }    // End fp
  }      // End f
}        // End function

void OpticalPumping_Method::apply_transverse_field_eg(DM_container *in) {
  // The energy difference between e & g states is so large that the
  // mixing is so small that this is definitely small enough to
  // ignore.
  in = in;
  /*
  for (int e = 0; e < numEStates; e++) {
    for (int g = 0; g < numGStates; g++) {
      // Get ready for transverse magnetic field!
      gsl_complex Bx = gsl_complex_rect(0.0, 0.0);
      gsl_complex left = gsl_complex_rect(0.0, 0.0);
      gsl_complex rigt = gsl_complex_rect(0.0, 0.0);
      // printf("For e = %d, comparing %d to %d\n", e, MFe2_Vector[e],
      //        MFe2_Vector[e+1]);
      // Is there a state with Mf` = Mf+1 and F` = F?
      if (e+1 < numEStates) {           // Keep things in bounds
        if (MFe2_Vector[e]+2 == MFe2_Vector[e+1]) {  // Sanity check
          gsl_complex tmp = in->eg[e+1][g];
          // printf("Bx = %g + %g i + ", GSL_REAL(tmp), GSL_IMAG(tmp));
          gsl_complex_mul_real(tmp, cPlus_E[e]);
          left = gsl_complex_add(left, tmp);
        }
      }
      // Is there a state with Mf` = Mf-1 and F` = F?
      if (e > 0) {                      // Keep things in bounds
        if (MFe2_Vector[e]-2 == MFe2_Vector[e-1]) {  // Sanity check
          gsl_complex tmp = in->eg[e-1][g];
          // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
          gsl_complex_mul_real(tmp, cMins_E[e]);
          left = gsl_complex_add(left, tmp);
        }
      }
      left = gsl_complex_mul_real(left, gFactor_E[e]);
      // Is there a state with Mf` = Mf+1 and F` = F?
      if (g+1 < numGStates) {           // Keep things in bounds
        if (MFg2_Vector[g]+2 == MFg2_Vector[g+1]) {  // Sanity check
          gsl_complex tmp = in->eg[e][g+1];
          // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
          gsl_complex_mul_real(tmp, cPlus_G[g]);
          rigt = gsl_complex_sub(rigt, tmp);
        }
      }
      // Is there a state with Mf` = Mf-1 and F` = F?
      if (g > 0) {                      // Keep things in bounds
        if (MFg2_Vector[g]-2 == MFg2_Vector[g-1]) {  // Sanity check
          gsl_complex tmp = in->eg[e][g-1];
          // printf("%g + %g i\n", GSL_REAL(tmp), GSL_IMAG(tmp));
          gsl_complex_mul_real(tmp, cMins_G[g]);
          rigt = gsl_complex_sub(rigt, tmp);
        }
      }
      rigt = gsl_complex_mul_real(rigt, gFactor_G);
      Bx = gsl_complex_add(left, rigt);
      Bx = gsl_complex_mul_imag(Bx,
                                -_bohr_magneton*eigen.field.B_x
                                /2.0/_planck_hbar);
      // printf("Bx contribution is %g + %g i\n", GSL_REAL(Bx), GSL_IMAG(Bx));
      dm_derivs->eg[e][g] = gsl_complex_add(dm_derivs->eg[e][g], Bx);
    }  // End g loop
  }    // End e loop
  */
}

void OpticalPumping_Method::apply_transverse_field_ef(DM_container *in) {
  // The energy difference between e & g states is so large that the
  // mixing is so small that this is definitely small enough to
  // ignore.
  in = in;
  /*
  for (int e = 0; e < numEStates; e++) {
    for (int f = 0; f < numFStates; f++) {
      // Get ready for transverse magnetic field!
      gsl_complex Bx = gsl_complex_rect(0.0, 0.0);
      gsl_complex left = gsl_complex_rect(0.0, 0.0);
      gsl_complex rigt = gsl_complex_rect(0.0, 0.0);
      // printf("For e = %d, comparing %d to %d\n", e, MFe2_Vector[e],
      //        MFe2_Vector[e+1]);
      // Is there a state with Mf` = Mf+1 and F` = F?
      if (e+1 < numEStates) {           // Keep things in bounds
        if (MFe2_Vector[e]+2 == MFe2_Vector[e+1]) {  // Sanity check
          gsl_complex tmp = in->ef[e+1][f];
          // printf("Bx = %g + %g i + ", GSL_REAL(tmp), GSL_IMAG(tmp));
          gsl_complex_mul_real(tmp, cPlus_E[e]);
          left = gsl_complex_add(left, tmp);
        }
      }
      // Is there a state with Mf` = Mf-1 and F` = F?
      if (e > 0) {                      // Keep things in bounds
        if (MFe2_Vector[e]-2 == MFe2_Vector[e-1]) {  // Sanity check
          gsl_complex tmp = in->ef[e-1][f];
          // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
          gsl_complex_mul_real(tmp, cMins_E[e]);
          left = gsl_complex_add(left, tmp);
        }
      }
      left = gsl_complex_mul_real(left, gFactor_E[e]);
      // Is there a state with Mf` = Mf+1 and F` = F?
      if (f+1 < numFStates) {           // Keep things in bounds
        if (MFf2_Vector[f]+2 == MFf2_Vector[f+1]) {  // Sanity check
          gsl_complex tmp = in->ef[e][f+1];
          // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
          gsl_complex_mul_real(tmp, cPlus_F[f]);
          rigt = gsl_complex_sub(rigt, tmp);
        }
      }
      // Is there a state with Mf` = Mf-1 and F` = F?
      if (f > 0) {                      // Keep things in bounds
        if (MFf2_Vector[f]-2 == MFf2_Vector[f-1]) {  // Sanity check
          gsl_complex tmp = in->ef[e][f-1];
          // printf("%g + %g i\n", GSL_REAL(tmp), GSL_IMAG(tmp));
          gsl_complex_mul_real(tmp, cMins_F[f]);
          rigt = gsl_complex_sub(rigt, tmp);
        }
      }
      rigt = gsl_complex_mul_real(rigt, gFactor_F);
      Bx = gsl_complex_add(left, rigt);
      Bx = gsl_complex_mul_imag(Bx,
                                -_bohr_magneton*eigen.field.B_x
                                /2.0/_planck_hbar);
      // printf("Bx contribution is %g + %g i\n", GSL_REAL(Bx), GSL_IMAG(Bx));
      dm_derivs->ef[e][f] = gsl_complex_add(dm_derivs->ef[e][f], Bx);
    }  // End f loop
  }    // End e loop
  */
}      // End function

void OpticalPumping_Method::apply_transverse_field_fg(DM_container *in) {
  // Working in the F-basis, each state has a definite F value and
  // therefore Bx cannot mix states with different F.
  in = in;
}      // End function

void OpticalPumping_Method::apply_transverse_field(DM_container *in) {
  apply_transverse_field_ee(in);
  apply_transverse_field_gg(in);
  apply_transverse_field_ff(in);
}


