// Copyright 2012 Benjamin Fenker
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_coupling.h>
#include "include/optical_pumping_method.h"
#include "include/units.h"
using std::vector;
OpticalPumping_Method::OpticalPumping_Method() {}
OpticalPumping_Method::OpticalPumping_Method(atom_data atom,
                                             magnetic_field_data field,
                                             Laser_data set_laser_fe,
                                             Laser_data set_laser_ge)
  : numEStates(atom.numEStates), numFStates(atom.numFStates),
    numGStates(atom.numGStates), tau(atom.tau), laser_ge(set_laser_ge),
    laser_fe(set_laser_fe), Fe2_Vector(numEStates, 0),
    MFe2_Vector(numEStates, 0), Ff2_Vector(numFStates, 0),
    MFf2_Vector(numFStates, 0), Fg2_Vector(numGStates, 0),
    MFg2_Vector(numGStates, 0), nu_E(numEStates, 0.0), nu_F(numFStates, 0.0),
    nu_G(numGStates, 0.0),
    a_eg(numEStates, vector<vector<double> >(numGStates,
                                             vector<double>(3, 0.0))),
    a_ef(numEStates, vector<vector<double> >(numFStates,
                                             vector<double>(3, 0.0))),
    rho_ee(numEStates, vector<gsl_complex>(numEStates,
                                           gsl_complex_rect(0.0, 0.0))),
    rho_ff(numFStates, vector<gsl_complex>(numFStates,
                                           gsl_complex_rect(0.0, 0.0))),
    rho_gg(numGStates, vector<gsl_complex>(numGStates,
                                           gsl_complex_rect(0.0, 0.0))),
    delta_ef(numEStates, vector<gsl_complex>(numFStates,
                                                 gsl_complex_rect(0.0, 0.0))),
    delta_eg(numEStates, vector<gsl_complex>(numGStates,
                                                 gsl_complex_rect(0.0, 0.0))),
    delta_fg(numFStates, vector<gsl_complex>(numGStates,
                                                 gsl_complex_rect(0.0, 0.0))) {
  printf("OpticalPumping_Method::OpticalPumping_Method(...)\n");

  for (int i = 0; i < 3; i++) {
    printf("%4.2G\t%4.2G\n", laser_ge.polarization[i],
           laser_fe.polarization[i]);
  }
  setup_quantum_numbers(atom);
  setup_frequencies_excited(atom, field);
  setup_frequencies_ground(atom, field);
  setup_eg_coupling(atom);
  setup_ef_coupling(atom);
  print_couplings(stdout);
  setup_pop_uniform_ground();
  // print_density_matrix(stdout);
  // print_data(stdout, 0.0);
  // for (int e = 0; e < numEStates; e++) printf("%15.10G\n",nu_E[e]);
}

OpticalPumping_Method::~OpticalPumping_Method() {}

void OpticalPumping_Method::setup_quantum_numbers(int I2, int Je2) {
  bool debug = true;
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

void OpticalPumping_Method::setup_frequencies_excited(int I2, int Je2,
                                                      double excitation,
                                                      double hyperfine_const,
                                                      double mu_B, double g_I,
                                                      double B_z) {
  for (int e = 0; e < numEStates; e++) {
    nu_E[e] = set_frequency(excitation, I2, Je2, Fe2_Vector[e], MFe2_Vector[e],
                            2, hyperfine_const, mu_B, g_I, B_z);
  }
}

void OpticalPumping_Method::setup_frequencies_excited(
                               atom_data atom,
                               magnetic_field_data field) {
  setup_frequencies_excited(atom.I2, atom.Je2, atom.nu_excited, atom.Aj_e,
                            field.mu_B, atom.g_I, field.B_z);
}

void OpticalPumping_Method::setup_frequencies_ground(int I2,
                                                     double hyperfine_const,
                                                     double mu_B, double g_I,
                                                     double B_z) {
  for (int g = 0; g < numGStates; g++) {
    nu_G[g] = set_frequency(0.0, I2, 1, Fg2_Vector[g], MFg2_Vector[g], 0,
                            hyperfine_const, mu_B, g_I, B_z);
  }
  for (int f = 0; f < numFStates; f++) {
    nu_F[f] = set_frequency(0.0, I2, 1, Ff2_Vector[f], MFf2_Vector[f], 0,
                            hyperfine_const, mu_B, g_I, B_z);
  }
}

void OpticalPumping_Method::setup_frequencies_ground(
                                atom_data atom,
                                magnetic_field_data field) {
  setup_frequencies_ground(atom.I2, atom.Aj_g, field.mu_B, atom.g_I, field.B_z);
}
double OpticalPumping_Method::set_frequency(double excitation, int I2, int J2,
                                            int F2, int Mf2, int L2,
                                            double hyperfine_const,
                                            double mu_B, double g_I,
                                            double B_z) {
  bool debug = true;
  Eigenvector_Helper alk;
  double g_f = alk.calc_gf(F2, J2, I2, L2, 1, g_I);

  // See DM equation 3.7
  double hyperfine = static_cast<double>((F2*(F2+2)) - (I2*(I2+2)) -
                                         (J2*(J2+2)));
  hyperfine *= hyperfine_const/8.0;
  if (debug) printf("Hyperfine: %15.10G MHz \t", hyperfine/_MHz);
  if (debug) printf("g_f = %6.4G, mu_B = %6.4G MHz/G, B_z = %6.4G G, ", g_f,
                    mu_B/(_MHz/_G), B_z/_G);
  double zeeman = static_cast<double>(Mf2);
  zeeman *= g_f * mu_B * B_z / 2.0;
  if (debug) printf("Zeeman: %15.10G MHz\t", zeeman/_MHz);

  if (debug) printf("Total: %15.10G MHz\n",
                    (excitation + hyperfine + zeeman)/_MHz);
  return (excitation + hyperfine + zeeman);
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

void OpticalPumping_Method::setup_pop_uniform_ground() {
  double startPop = numFStates + numGStates;
  startPop = 1.0 / startPop;
  for (int f = 0; f < numFStates; f++) GSL_SET_REAL(&rho_ff[f][f], startPop);
  for (int g = 0; g < numGStates; g++) GSL_SET_REAL(&rho_gg[g][g], startPop);
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
  for (int g = 0; g < numGStates; g++) sum += GSL_REAL(rho_gg[g][g]);
  for (int f = 0; f < numFStates; f++) sum += GSL_REAL(rho_ff[f][f]);
  for (int e = 0; e < numEStates; e++) sum += GSL_REAL(rho_ee[e][e]);
  return sum;
}

double OpticalPumping_Method::get_polarization() { return 0.0; }

double OpticalPumping_Method::get_alignment() { return 0.0; }

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
    if (fabs(GSL_IMAG(rho_gg[g][g])) > eps) hermit = false;
  }
  for (int f = 0; f < numFStates; f++) {
    if (fabs(GSL_IMAG(rho_ff[f][f])) > eps) hermit = false;
  }
  for (int e = 0; e < numEStates; e++) {
    if (fabs(GSL_IMAG(rho_ee[e][e])) > eps) hermit = false;
  }

  // Now check that the off-diagonals have the appropriate symmetry
  for (int g = 0; g < numGStates; g++) {
    for (int gp = 0; gp < numGStates; gp++) {
      double imaginary_sum = fabs(GSL_IMAG(rho_gg[g][gp]) +
                                   GSL_IMAG(rho_gg[gp][g]));
      double real_diff = fabs(GSL_REAL(rho_gg[g][gp]) -
                               GSL_REAL(rho_gg[gp][g]));
      if ((imaginary_sum > eps) || (real_diff > eps)) hermit = false;
    }
  }
  for (int f = 0; f < numFStates; f++) {
    for (int fp = 0; fp < numFStates; fp++) {
      double imaginary_sum = fabs(GSL_IMAG(rho_ff[f][fp]) +
                                   GSL_IMAG(rho_ff[fp][f]));
      double real_diff = fabs(GSL_REAL(rho_ff[f][fp]) -
                               GSL_REAL(rho_ff[fp][f]));
      if ((imaginary_sum > eps) || (real_diff > eps)) hermit = false;
    }
  }
  for (int e = 0; e < numEStates; e++) {
    for (int ep = 0; ep < numGStates; ep++) {
      double imaginary_sum = fabs(GSL_IMAG(rho_ee[e][ep]) +
                                   GSL_IMAG(rho_ee[ep][e]));
      double real_diff = fabs(GSL_REAL(rho_ee[e][ep]) -
                               GSL_REAL(rho_ee[ep][e]));
      if ((imaginary_sum > eps) || (real_diff > eps)) hermit = false;
    }
  }
  return hermit;
}

void OpticalPumping_Method::print_data(FILE *des, double time) {
  double pop = get_total_population();
  double pol = get_polarization();
  double ali = get_alignment();

  fprintf(des, "%8.6G\t", time/_us);
  for (int g = 0; g < numGStates; g++) fprintf(des, "%8.6G\t",
                                               GSL_REAL(rho_gg[g][g]));
  for (int f = 0; f < numFStates; f++) fprintf(des, "%8.6G\t",
                                               GSL_REAL(rho_ff[f][f]));
  for (int e = 0; e < numEStates; e++) fprintf(des, "%8.6G\t",
                                               GSL_REAL(rho_ee[e][e]));
  fprintf(des, "%8.6G\t%8.6G\t%8.6G", pop, pol, ali);
  fprintf(des, "\n");
}

void OpticalPumping_Method::print_density_matrix(FILE *des) {
  // This function is useful when printing out the whole density matrix at one
  // time step to the screen.  It prints the real part as one matrix, tabs over
  // and prints the imaginary part as another matrix to the right of them first}
  fprintf(des, "rho_ee\n");
  for (int r = 0; r < numEStates; r++) {
    for (int c = 0; c < numEStates; c++) {
      fprintf(des, "%10.6G   ", GSL_REAL(rho_ee[r][c]));
    }
    fprintf(des, "\t\t");
    for (int c = 0; c < numEStates; c++) {
      fprintf(des, "%10.6G   ", GSL_IMAG(rho_ee[r][c]));
    }
    fprintf(des, "\n");
  }
  fprintf(des, "rho_gg\n");
  for (int r = 0; r < numGStates; r++) {
    for (int c = 0; c < numGStates; c++) {
      fprintf(des, "%10.6G   ", GSL_REAL(rho_gg[r][c]));
    }
    fprintf(des, "\t\t");
    for (int c = 0; c < numGStates; c++) {
      fprintf(des, "%10.6G   ", GSL_IMAG(rho_gg[r][c]));
    }
    fprintf(des, "\n");
  }
  fprintf(des, "rho_ff\n");
  for (int r = 0; r < numFStates; r++) {
    for (int c = 0; c < numFStates; c++) {
      fprintf(des, "%10.6G   ", GSL_REAL(rho_ff[r][c]));
    }
    fprintf(des, "\t\t");
    for (int c = 0; c < numFStates; c++) {
      fprintf(des, "%10.6G   ", GSL_IMAG(rho_ff[r][c]));
    }
    fprintf(des, "\n");
  }
  fprintf(des, "delta_rho_eg\n");
  for (int r = 0; r < numEStates; r++) {
    for (int c = 0; c < numGStates; c++) {
      fprintf(des, "%10.6G   ", GSL_REAL(delta_eg[r][c]));
    }
    fprintf(des, "\t\t");
    for (int c = 0; c < numGStates; c++) {
      fprintf(des, "%10.6G   ", GSL_IMAG(delta_eg[r][c]));
    }
    fprintf(des, "\n");
  }
  fprintf(des, "delta_rho_ef\n");
  for (int r = 0; r < numEStates; r++) {
    for (int c = 0; c < numFStates; c++) {
      fprintf(des, "%10.6G   ", GSL_REAL(delta_ef[r][c]));
    }
    fprintf(des, "\t\t");
    for (int c = 0; c < numFStates; c++) {
      fprintf(des, "%10.6G   ", GSL_IMAG(delta_ef[r][c]));
    }
    fprintf(des, "\n");
  }
  fprintf(des, "delta_fg\n");
  for (int r = 0; r < numFStates; r++) {
    for (int c = 0; c < numGStates; c++) {
      fprintf(des, "%10.6G   ", GSL_REAL(delta_fg[r][c]));
    }
    fprintf(des, "\t\t");
    for (int c = 0; c < numGStates; c++) {
      fprintf(des, "%10.6G   ", GSL_IMAG(delta_fg[r][c]));
    }
    fprintf(des, "\n");
  }
}

