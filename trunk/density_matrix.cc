// Copyright 2012 Benjamin Fenker

#include <stdio.h>
#include <math.h>
#include <vector>
#include <omp.h>
#include "include/density_matrix.h"
#include "include/units.h"

using std::vector;
using std::min;
Density_Matrix::Density_Matrix(atom_data atom, magnetic_field_data field,
                               Laser_data set_laser_fe, Laser_data set_laser_ge,
                               coherence_flags flags)
  : OpticalPumping_Method(atom, field, set_laser_fe, set_laser_ge),
    dipole_moment_eg(numEStates, vector<double>(numGStates, 0.0)),
    dipole_moment_ef(numEStates, vector<double>(numFStates, 0.0)),
    drho_ee(numEStates, vector<gsl_complex>(numEStates,
                                            gsl_complex_rect(0.0, 0.0))),
    drho_ff(numFStates, vector<gsl_complex>(numFStates,
                                            gsl_complex_rect(0.0, 0.0))),
    drho_gg(numGStates, vector<gsl_complex>(numGStates,
                                            gsl_complex_rect(0.0, 0.0))),
    ddelta_ef(numEStates, vector<gsl_complex>(numFStates,
                                            gsl_complex_rect(0.0, 0.0))),
    ddelta_eg(numEStates, vector<gsl_complex>(numGStates,
                                            gsl_complex_rect(0.0, 0.0))),
    ddelta_fg(numFStates, vector<gsl_complex>(numGStates,
                                            gsl_complex_rect(0.0, 0.0))) {
  printf("DensityMatrix::Density_Matrix(...)\n");
  es_Zeeman = flags.zCoherences;
  gs_Zeeman = flags.zCoherences;
  es_hyperfine = flags.hfCoherences_ex;
  setup_dipole_moments(atom.linewidth);
  print_rabi_frequencies(stdout);

  totalTerms = (atom.numFStates + atom.numGStates + atom.numEStates);
  totalTerms = 2*totalTerms*totalTerms;
  data.gStart = 0;
  data.gEnd = data.numGStates * data.numGStates * 2;
  data.fStart = data.gEnd;
  data.fEnd = data.fStart + (data.numFStates * data.numFStates * 2);
  data.eStart = data.fEnd;
  data.eEnd = data.eStart + (data.numEStates * data.numEStates * 2);
  data.dipole_moment_eg = dipole_moment_eg;
  data.dipole_moment_ef = dipole_moment_ef;

  population = new double[totalTerms];
  // This is the ordering of the density matrix terms:

  // First is the rho_gg portion of the matrix going across in one row before
  // moving down to the second column, etc.  Real then imaginary term by term
  // For example population[0] = real(rho_gg[0][0]) and
  // population[1] = imag(rho_gg[0][1]) and population[2] = real(rho_gg[0][1])

  // After this the same procedure is repeated continuing with rho_ff
  // immediately following the last entry in the rho_gg matrix.  Immediately
  // follwing this is the rho_ee matrix.

  // Next is the coherences, starting with ef then eg then fg.  The ordering
  // for these is the same as above:  real then imaginary for each term
  // advancing across a row before moving to the next column
  double totalStates = data.numFStates + data.numGStates + data.numEStates;
  for (int i = data.gStart; i < data.eEnd; i++) population[i] = 0.0;
  for (int g = data.gStart; g < data.gEnd; g++) {
      if (g % ((data.atom.numGStates+1)*2) == 0) {
        printf("setting population[%d] = 1\n", g);
        population[g] = 1.0/totalStates;
      } else {
        population[g] = 0.0;
    }
  }
  for (int f = data.fStart; f < data.fEnd; f++) {
    if ((f-data.fStart) % ((data.atom.numFStates+1)*2) == 0) {
      printf("setting population[%d] = 1\n", f);
      population[f] = 1.0/totalStates;
    } else {
      population[f] = 0.0;
    }
  }
  for (int e = data.eStart; e < data.eEnd; e++) {
    if ((e-data.eStart) % ((data.atom.numEStates+1)*2) == 0) {
      printf("setting population[%d] = 1\n", e);
      population[e] = 0.0;
    } else {
      population[e] = 0.0;
    }
  }
  for (int i = data.eEnd; i < totalTerms; i++) population[i] = 0.0;
}

Density_Matrix::~Density_Matrix() {
  delete population;
}

int Density_Matrix::update_population_gsl(double t, const double y[],
                                           double deriv[], void *params) {
  //  printf("Density_Matrix::update_population_gsl(...)\n");
  op_data_for_gsl data = *(reinterpret_cast<op_data_for_gsl *>(params));
  //  deriv[0] = -y[0] / data.atom.tau;

  // First re-pack the y array into gsl_complexes
  // Also make and initalize to zero matrices of gsl_complex to hold the
  // derivatives.  The point of all this is to copy and paste the physics code
  // I've already built and to let the computer handle the complex arithmatic
  gsl_complex **rho_gg = new gsl_complex*[data.numGStates];
  gsl_complex **rho_ff = new gsl_complex*[data.numFStates];
  gsl_complex **rho_ee = new gsl_complex*[data.numEStates];
  gsl_complex **rho_ef = new gsl_complex*[data.numEStates];
  gsl_complex **rho_eg = new gsl_complex*[data.numEStates];
  gsl_complex **rho_fg = new gsl_complex*[data.numFStates];

  gsl_complex **drho_gg = new gsl_complex*[data.numGStates];
  gsl_complex **drho_ff = new gsl_complex*[data.numFStates];
  gsl_complex **drho_ee = new gsl_complex*[data.numEStates];
  gsl_complex **drho_ef = new gsl_complex*[data.numEStates];
  gsl_complex **drho_eg = new gsl_complex*[data.numEStates];
  gsl_complex **drho_fg = new gsl_complex*[data.numFStates];

  int index, extra;
  for (int i = 0; i < data.numEStates; i++) {
    rho_ee[i] = new gsl_complex[data.numEStates];
    drho_ee[i] = new gsl_complex[data.numEStates];
    rho_ef[i] = new gsl_complex[data.numFStates];
    drho_ef[i] = new gsl_complex[data.numFStates];
    rho_eg[i] = new gsl_complex[data.numGStates];
    drho_eg[i] = new gsl_complex[data.numGStates];
    if (i < data.numFStates) {
      rho_ff[i] = new gsl_complex[data.numFStates];
      drho_ff[i] = new gsl_complex[data.numFStates];
      rho_fg[i] = new gsl_complex[data.numGStates];
      drho_fg[i] = new gsl_complex[data.numGStates];
    }
    if (i < data.numGStates) {
      rho_gg[i] = new gsl_complex[data.numGStates];
      drho_gg[i] = new gsl_complex[data.numGStates];
    }
    for (int j = 0; j < data.numEStates; j++) {
      index = (2*((i*data.numEStates) + j)) + data.eStart;
      rho_ee[i][j] = gsl_complex_rect(y[index], y[index+1]);
      drho_ee[i][j] = gsl_complex_rect(0.0, 0.0);
      // printf("rho_ee[%d][%d] = %5.3G + %5.3Gi\n", i, j,
      // GSL_REAL(rho_ee[i][j]), GSL_IMAG(rho_ee[i][j]));
      if (j < data.numFStates) {
        index = (2*((i*data.numFStates) + j)) + data.eEnd;
        rho_ef[i][j] = gsl_complex_rect(y[index], y[index+1]);
        drho_ef[i][j] = gsl_complex_rect(0.0, 0.0);
        // printf("rho_ef[%d][%d] = %5.3G + %5.3Gi\n", i, j,
        //       GSL_REAL(rho_ef[i][j]), GSL_IMAG(rho_ef[i][j]));
        if (i < data.numFStates && j < data.numGStates) {
          extra = 2*(data.numEStates * data.numFStates);
          extra += 2*(data.numEStates * data.numGStates);
          index = (2*((i*data.numGStates) + j)) + data.eEnd + extra;
          rho_fg[i][j] = gsl_complex_rect(y[index], y[index+1]);
          drho_fg[i][j] = gsl_complex_rect(0.0, 0.0);
          // printf("rho_fg[%d][%d] = %5.3G + %5.3Gi\n", i, j,
          //       GSL_REAL(rho_fg[i][j]), GSL_IMAG(rho_fg[i][j]));
        }
        if (i < data.numFStates && j < data.numFStates) {
          index = (2*((i*data.numFStates) + j)) + data.fStart;
          rho_ff[i][j] = gsl_complex_rect(y[index], y[index+1]);
          drho_ff[i][j] = gsl_complex_rect(0.0, 0.0);
          // printf("rho_ff[%d][%d] = %5.3G + %5.3Gi\n", i, j,
          //     GSL_REAL(rho_ff[i][j]), GSL_IMAG(rho_ff[i][j]));
        }
      }
      if (j < data.numGStates) {
        extra = data.numEStates * data.numFStates * 2;
        index = (2*((i*data.numGStates) + j)) + data.eEnd + extra;
        rho_eg[i][j] = gsl_complex_rect(y[index], y[index+1]);
        drho_eg[i][j] = gsl_complex_rect(0.0, 0.0);
        // printf("rho_eg[%d][%d] = %5.3G + %5.3Gi\n", i, j,
        //   GSL_REAL(rho_eg[i][j]), GSL_IMAG(rho_eg[i][j]));
        if (i < data.numGStates && j < data.numGStates) {
          index = (2*((i*data.numGStates) + j)) + data.gStart;
          rho_gg[i][j] = gsl_complex_rect(y[index], y[index+1]);
        drho_gg[i][j] = gsl_complex_rect(0.0, 0.0);
        // printf("rho_gg[%d][%d] = %5.3G + %5.3Gi\n", i, j,
        //     GSL_REAL(rho_gg[i][j]), GSL_IMAG(rho_gg[i][j]));
        }
      }
    }
  }
  // Finally done extracting the easy-to-work with GSL_COMPLEX from the
  // dmenaded-by-gsl doubles.  Now time to do some physics!

  // The bare minimum is to include only diagonal matrix elements as well
  // as the optical coherences.  Spontaneous decay will be included as well

  // Booleans that eventually I need to have passed from on high
  bool zCoherences = false, hfCoherences_ex = false, hfCoherences_gr = false;
  bool es_Zeeman = false, gs_Zeeman;
  // Excited states rho_ee:  Equation 32
  // Zeeman coherences are off-diagonals where F = F`
  // Hyperfine coherences are off-diagonal where M_f = M_f`
  
  for (int e = 0; e < data.numEStates; e++) {
    for (int ep = 0; ep < data.numEStates; ep++) {
      // if ((e == ep) || (es_Zeeman && MFe2_Vector[e] == MFe2_Vector[ep]) ||
      // (es_hyperfine && Fe2_Vector[e] == Fe2_Vector[ep]) ||
      // (es_hyperfine && es_Zeeman)) {
      if (e == ep) {            // Need a way to flag coherences on and off
        { // G-laser term
          gsl_complex q_sum = gsl_complex_rect(0.0, 0.0);
          for (int q = 0; q < 3; q++) {
            gsl_complex g_sum = gsl_complex_rect(0.0, 0.0);
            for (int gpp = 0; gpp < data.numGStates; gpp++) {
              gsl_complex left = gsl_complex_rect(data.dipole_moment_eg[e][gpp],
                                                  0.0);
              left = gsl_complex_mul_real(left, data.a_eg[e][gpp][q]);
              left = gsl_complex_mul(left,
                                     gsl_complex_conjugate(rho_eg[ep][gpp]));
              gsl_complex right = gsl_complex_rect(pow(-1.0, (0))*
                                      data.dipole_moment_eg[ep][gpp],
                                                   0.0);
              right = gsl_complex_mul_real(right, data.a_eg[ep][gpp][q]);
              right = gsl_complex_mul(right,
                                      rho_eg[e][gpp]);
              left = gsl_complex_sub(left, right);
              g_sum = gsl_complex_add(g_sum, left);
            }
            g_sum = gsl_complex_mul_real(g_sum, data.laser_ge.polarization[q]);
            q_sum = gsl_complex_add(q_sum, g_sum);
          }
          q_sum = gsl_complex_mul_real(q_sum, data.laser_ge.field);
          q_sum = gsl_complex_mul_imag(q_sum, 0.5/ _planck_hbar);
          drho_ee[e][ep] = gsl_complex_add(drho_ee[e][ep], q_sum);
        }
      }
    }
  }
  
  // G-Ground states rho_gg:  Equation 34
  // Zeeman coherences are when g != g`
  for (int g = 0; g < data.numGStates; g++) {
    for (int gp = 0; gp < data.numGStates; gp++) {
      gsl_complex oCoherences = gsl_complex_rect(0.0, 0.0);
      if ((g == gp) || gs_Zeeman) {
        for (int q = 0; q < 3; q++) {
          gsl_complex qTerm = gsl_complex_rect(0.0, 0.0);
          for (int epp = 0; epp < data.numEStates; epp++) {
            gsl_complex left, right;  // First and second halves in
            // top line of equation 34
            double dTemp = pow(-1.0, 0) * data.dipole_moment_eg[epp][g];
            dTemp *= data.a_eg[epp][g][gp];
            left = gsl_complex_rect(dTemp, 0.0);
            left = gsl_complex_mul(left,
                                   rho_eg[epp][gp]);
            dTemp = data.dipole_moment_eg[epp][gp] * data.a_eg[epp][gp][q];
            right = gsl_complex_rect(dTemp, 0.0);
            right = gsl_complex_mul(right,
                                    gsl_complex_conjugate(rho_eg[epp][g]));
            left = gsl_complex_sub(left, right);
            qTerm = gsl_complex_add(qTerm, left);
          }
          qTerm = gsl_complex_mul_real(qTerm, data.laser_ge.polarization[q]);
          oCoherences = gsl_complex_add(oCoherences, qTerm);
        }
        double mulTerm = data.laser_ge.field / (2*_planck_hbar);
        oCoherences = gsl_complex_mul_imag(oCoherences, mulTerm);
        drho_gg[g][gp] = gsl_complex_add(drho_gg[g][gp], oCoherences);
      }
    }
  }

  // EG-Optical Coherences rho_eg:  Equation 36
  bool debug_eg = false;
  for (int e = 0; e < data.numEStates; e++) {
    for (int g = 0; g < data.numGStates; g++) {
      if (debug_eg) printf("e = %d, g = %d\n", e, g);
      gsl_complex gLaser_term = gsl_complex_rect(0.0, 0.0);
      for (int q = 0; q < 3; q++) {
        if (data.laser_ge.polarization[q] > 0.0) {
          if (debug_eg) printf("\t q = %d\n", q);
          gsl_complex gsum = gsl_complex_rect(0.0, 0.0);
          gsl_complex esum = gsl_complex_rect(0.0, 0.0);
          for (int gpp = 0; gpp < data.numGStates; gpp++) {
            gsl_complex temp = gsl_complex_rect(data.dipole_moment_eg[e][gpp],
                                                0.0);
            temp = gsl_complex_mul_real(temp, data.a_eg[e][gpp][q]);
            temp = gsl_complex_mul(temp, rho_gg[g][gpp]);
            gsum = gsl_complex_add(gsum, temp);
          }
          for (int epp = 0; epp < data.numEStates; epp++) {
            gsl_complex temp = gsl_complex_rect(data.dipole_moment_eg[epp][g],
                                                0.0);
            temp = gsl_complex_mul_real(temp, data.a_eg[epp][g][q]);
            temp = gsl_complex_mul(temp, rho_ee[e][epp]);
            esum = gsl_complex_add(esum, temp);
          }
          // printf("gsum = %8.6G + %8.6Gi\n", GSL_REAL(gsum), GSL_IMAG(gsum));
          // printf("esum = %8.6G + %8.6Gi\n", GSL_REAL(esum), GSL_IMAG(esum));
          gsum = gsl_complex_sub(gsum, esum);
          gsum = gsl_complex_mul_real(gsum, data.laser_ge.polarization[q]);
          gLaser_term = gsl_complex_add(gLaser_term, gsum);
        }
      }
      gLaser_term = gsl_complex_mul_real(gLaser_term, data.laser_ge.field);
      gLaser_term = gsl_complex_mul_imag(gLaser_term, 1.0/(2.0*_planck_hbar));
      drho_eg[e][g] = gsl_complex_add(drho_eg[e][g], gLaser_term);

      // Detuning term -i(\omega_eg - \omega_1)*\rho_eg
      double detune_d = ((data.nu_E[e] - data.nu_G[g]) - data.laser_ge.nu);
      detune_d *= -2.0*M_PI;
      gsl_complex detune_g = gsl_complex_rect(0.0, detune_d);
      detune_g = gsl_complex_mul(detune_g, rho_eg[e][g]);
      drho_eg[e][g] = gsl_complex_add(drho_eg[e][g], detune_g);
    }
  }
  
  // Print out the density matrix
  
  // After working with the gsl complex, repack the drhos into the derivatives
  // that GSL_ODE needs
  index = 0;
  for (int g = 0; g < data.numGStates; g++) {
    for (int gp = 0; gp < data.numGStates; gp++) {
      deriv[index] = GSL_REAL(drho_gg[g][gp]);
      deriv[index+1] = GSL_IMAG(drho_gg[g][gp]);
      index += 2;
    }
  }
  for (int f  = 0; f < data.numFStates; f++) {
    for (int fp = 0; fp < data.numFStates; fp++) {
      deriv[index] = GSL_REAL(drho_ff[f][fp]);
      deriv[index+1] = GSL_IMAG(drho_ff[f][fp]);
      index += 2;
    }
  }
  for (int e  = 0; e < data.numEStates; e++) {
    for (int ep = 0; ep < data.numEStates; ep++) {
      deriv[index] = GSL_REAL(drho_ee[e][ep]);
      deriv[index+1] = GSL_IMAG(drho_ee[e][ep]);
      index += 2;
    }
  }
  for (int e = 0; e < data.numEStates; e++) {
    for (int f = 0; f < data.numFStates; f++) {
      deriv[index] = GSL_REAL(drho_ef[e][f]);
      deriv[index+1] = GSL_IMAG(drho_ef[e][f]);
      index += 2;
    }
  }
  for (int e = 0; e < data.numEStates; e++) {
    for (int g = 0; g < data.numGStates; g++) {
      //      printf("drho_eg[%d][%d] = %8.6G + %8.6Gi\n", e, g,
      //     GSL_REAL(drho_eg[e][g]), GSL_IMAG(drho_eg[e][g]));
      deriv[index] = GSL_REAL(drho_eg[e][g]);
      deriv[index+1] = GSL_IMAG(drho_eg[e][g]);
      index += 2;
    }
  }
  for (int f = 0; f < data.numFStates; f++) {
    for (int g = 0; g < data.numGStates; g++) {
      deriv[index] = GSL_REAL(drho_fg[f][g]);
      deriv[index+1] = GSL_IMAG(drho_fg[f][g]);
      index += 2;
    }
  }
  //  printf("Finished on: %d\n", index);

  // Delete all the dynamically allocated arrays!
  for (int i = 0; i < data.numEStates; i++) {
    delete[] rho_ee[i];
    delete[] drho_ee[i];
    delete[] rho_ef[i];
    delete[] drho_ef[i];
    delete[] rho_eg[i];
    delete[] drho_eg[i];
    if (i < data.numGStates) {
      delete[] rho_gg[i];
      delete[] drho_gg[i];
    }
    if (i < data.numFStates) {
      delete[] rho_ff[i];
      delete[] drho_ff[i];
      delete[] rho_fg[i];
      delete[] drho_fg[i];
    }
  }
  return GSL_SUCCESS;
}

void Density_Matrix::update_population(double dt) {
  reset_dPop();
  // The bare minimum is to include only diagonal matrix elements as well
  // as the optical coherences.  Spontaneous decay will be included as well
  /*
  int e_check = 3;
  int g_check = 0;
  int q_check = 2;
  printf("D_eg[%d][%d] = %8.6G e*nm\t", e_check, g_check,
         dipole_moment_eg[e_check][g_check]/(_elementary_charge*_nm));
  printf("a_eg[%d][%d] = %8.6G\t", e_check, g_check,
         a_eg[e_check][g_check][q_check]);
  printf("g->e field = %8.6G V/m\t", laser_ge.field/(_V/_m));
  printf("hbar = %8.6G Js\n", _planck_hbar/(_J*_s));
  */
  // Excited states rho_ee:  Equation 32
  // Zeeman coherences are off-diagonals where F = F`
  // Hyperfine coherences are off-diagonal where M_f = M_f`
#pragma omp parallel 
  {
  
  for (int e = 0; e < numEStates; e++) {
    for (int ep = 0; ep < numEStates; ep++) {
      if ((e == ep) || (es_Zeeman && MFe2_Vector[e] == MFe2_Vector[ep]) ||
          (es_hyperfine && Fe2_Vector[e] == Fe2_Vector[ep]) ||
          (es_hyperfine && es_Zeeman)) {
        { // G-laser term
          gsl_complex q_sum = gsl_complex_rect(0.0, 0.0);
          for (int q = 0; q < 3; q++) {
            gsl_complex g_sum = gsl_complex_rect(0.0, 0.0);
            for (int gpp = 0; gpp < numGStates; gpp++) {
              gsl_complex left = gsl_complex_rect(dipole_moment_eg[e][gpp],
                                                  0.0);
              left = gsl_complex_mul_real(left, a_eg[e][gpp][q]);
              left = gsl_complex_mul(left,
                                     gsl_complex_conjugate(delta_eg[ep][gpp]));
              gsl_complex right = gsl_complex_rect(pow(-1.0, (0))*
                                                   dipole_moment_eg[ep][gpp],
                                                   0.0);
              right = gsl_complex_mul_real(right, a_eg[ep][gpp][q]);
              right = gsl_complex_mul(right,
                                      delta_eg[e][gpp]);
              left = gsl_complex_sub(left, right);
              g_sum = gsl_complex_add(g_sum, left);
            }
            g_sum = gsl_complex_mul_real(g_sum, laser_ge.polarization[q]);
            q_sum = gsl_complex_add(q_sum, g_sum);
          }
          q_sum = gsl_complex_mul_real(q_sum, laser_ge.field);
          q_sum = gsl_complex_mul_imag(q_sum, 0.5/ _planck_hbar);
          drho_ee[e][ep] = gsl_complex_add(drho_ee[e][ep], q_sum);
        } // End G-Laser term (no loop or if, just an arbitrary block
        { // F-laser term
          gsl_complex q_sum = gsl_complex_rect(0.0, 0.0);
          for (int q = 0; q < 3; q++) {
            gsl_complex f_sum = gsl_complex_rect(0.0, 0.0);
            for (int fpp = 0; fpp < numFStates; fpp++) {
              gsl_complex left = gsl_complex_rect(dipole_moment_ef[e][fpp],
                                                  0.0);
              left = gsl_complex_mul_real(left, a_ef[e][fpp][q]);
              left = gsl_complex_mul(left,
                                     gsl_complex_conjugate(delta_ef[ep][fpp]));
              gsl_complex right = gsl_complex_rect(pow(-1.0, (0))*
                                                   dipole_moment_ef[ep][fpp],
                                                   0.0);
              right = gsl_complex_mul_real(right, a_ef[ep][fpp][q]);
              right = gsl_complex_mul(right,
                                      delta_ef[e][fpp]);
              left = gsl_complex_sub(left, right);
              f_sum = gsl_complex_add(f_sum, left);
            }
            f_sum = gsl_complex_mul_real(f_sum, laser_fe.polarization[q]);
            q_sum = gsl_complex_add(q_sum, f_sum);
          }
          q_sum = gsl_complex_mul_real(q_sum, laser_fe.field);
          q_sum = gsl_complex_mul_imag(q_sum, 0.5/ _planck_hbar);
          drho_ee[e][ep] = gsl_complex_add(drho_ee[e][ep], q_sum);
        } // End F-Laser term (no loop or if, just an arbitrary block

      }   // End the if(coherences)
    }     // End ep loop
    // END THE EE PORTION OF THE DENSITY MATRIX

  // EG-Optical Coherences delta_eg:  Equation 36
    bool debug_eg = false;
    for (int g = 0; g < numGStates; g++) {
      if (debug_eg) printf("e = %d, g = %d\n", e, g);
      gsl_complex gLaser_term = gsl_complex_rect(0.0, 0.0);
      for (int q = 0; q < 3; q++) {
        if (laser_ge.polarization[q] > 0.0) {
          if (debug_eg) printf("\t q = %d\n", q);
          gsl_complex gsum = gsl_complex_rect(0.0, 0.0);
          gsl_complex esum = gsl_complex_rect(0.0, 0.0);
          for (int gpp = 0; gpp < numGStates; gpp++) {
            gsl_complex temp = gsl_complex_rect(dipole_moment_eg[e][gpp], 0.0);
            temp = gsl_complex_mul_real(temp, a_eg[e][gpp][q]);
            temp = gsl_complex_mul(temp, rho_gg[gpp][g]);
            gsum = gsl_complex_add(gsum, temp);
          }
          for (int epp = 0; epp < numEStates; epp++) {
            gsl_complex temp = gsl_complex_rect(dipole_moment_eg[epp][g], 0.0);
            temp = gsl_complex_mul_real(temp, a_eg[epp][g][q]);
            temp = gsl_complex_mul(temp, rho_ee[e][epp]);
            esum = gsl_complex_add(esum, temp);
          }
          gsum = gsl_complex_sub(gsum, esum);
          gsum = gsl_complex_mul_real(gsum, laser_ge.polarization[q]);
          gLaser_term = gsl_complex_add(gLaser_term, gsum);
        }
      }
      gLaser_term = gsl_complex_mul_real(gLaser_term, laser_ge.field);
      gLaser_term = gsl_complex_mul_imag(gLaser_term, 1.0/(2.0*_planck_hbar));
      ddelta_eg[e][g] = gsl_complex_add(ddelta_eg[e][g], gLaser_term);
      
      // Detuning term -i(\omega_eg - \omega_1)*\delta_eg
      double detune_d = -2*M_PI*((nu_E[e] - nu_G[g]) - laser_ge.nu);
      gsl_complex detune_g = gsl_complex_rect(0.0, detune_d);
      detune_g = gsl_complex_mul(detune_g, delta_eg[e][g]);
      ddelta_eg[e][g] = gsl_complex_add(ddelta_eg[e][g], detune_g);
    } // End g loop for eg-coherences

    // EF-Optical Coherences delta_ef:  Equation 35
    bool debug_ef = false;
    for (int f = 0; f < numFStates; f++) {
      if (debug_ef) printf("e = %d, f = %d\n", e, f);
      gsl_complex fLaser_term = gsl_complex_rect(0.0, 0.0);
      for (int q = 0; q < 3; q++) {
        if (laser_fe.polarization[q] > 0.0) {
          if (debug_ef) printf("\t q = %d\n", q);
          gsl_complex fsum = gsl_complex_rect(0.0, 0.0);
          gsl_complex esum = gsl_complex_rect(0.0, 0.0);
          for (int fpp = 0; fpp < numFStates; fpp++) {
            gsl_complex temp = gsl_complex_rect(dipole_moment_ef[e][fpp], 0.0);
            temp = gsl_complex_mul_real(temp, a_ef[e][fpp][q]);
            temp = gsl_complex_mul(temp, rho_ff[fpp][f]);
            fsum = gsl_complex_add(fsum, temp);
          }
          if (debug_ef) printf("fsum = %8.6G + %8.6Gi\n", GSL_REAL(fsum),
                               GSL_IMAG(fsum));
          for (int epp = 0; epp < numEStates; epp++) {
            gsl_complex temp = gsl_complex_rect(dipole_moment_ef[epp][f], 0.0);
            temp = gsl_complex_mul_real(temp, a_ef[epp][f][q]);
            temp = gsl_complex_mul(temp, rho_ee[e][epp]);
            esum = gsl_complex_add(esum, temp);
          }
          if (debug_ef) printf("esum = %8.6G + %8.6Gi\n", GSL_REAL(esum),
                               GSL_IMAG(esum));
          fsum = gsl_complex_sub(fsum, esum);
          fsum = gsl_complex_mul_real(fsum, laser_fe.polarization[q]);
          fLaser_term = gsl_complex_add(fLaser_term, fsum);
        }
      }
      fLaser_term = gsl_complex_mul_real(fLaser_term, laser_fe.field);
      fLaser_term = gsl_complex_mul_imag(fLaser_term, 1.0/(2.0*_planck_hbar));
      ddelta_ef[e][f] = gsl_complex_add(ddelta_ef[e][f], fLaser_term);

      // Detuning term -i(\omega_ef - \omega_1)*\delta_ef
      double detune_d = -2*M_PI*((nu_E[e] - nu_F[f]) - laser_fe.nu);
      gsl_complex detune_f = gsl_complex_rect(0.0, detune_d);
      detune_f = gsl_complex_mul(detune_f, delta_ef[e][f]);
      ddelta_ef[e][f] = gsl_complex_add(ddelta_ef[e][f], detune_f);
    } // End f loop for ef-coherences
  }       // End e loop
  // F-Ground states rho_ff:  Equation 34
  // Zeeman coherences are when f != f`
  for (int f = 0; f < numFStates; f++) {
    for (int fp = 0; fp < numFStates; fp++) {
      gsl_complex oCoherences = gsl_complex_rect(0.0, 0.0);
      if ((f == fp) || gs_Zeeman) {
        for (int q = 0; q < 3; q++) {
          gsl_complex qTerm = gsl_complex_rect(0.0, 0.0);
          for (int epp = 0; epp < numEStates; epp++) {
            gsl_complex left, right;  // First and second halves in
            // top line of equation 34
            left = gsl_complex_rect(pow(-1.0, (0))*
                                    dipole_moment_ef[epp][f]*a_ef[epp][f][q],
                                    0.0);
            left = gsl_complex_mul(left,
                                   delta_ef[epp][fp]);
            right = gsl_complex_rect(dipole_moment_ef[epp][fp]*a_ef[epp][fp][q],
                                     0.0);
            right = gsl_complex_mul(right,
                                    gsl_complex_conjugate(delta_ef[epp][f]));
            left = gsl_complex_sub(left, right);
            qTerm = gsl_complex_add(qTerm, left);
          }
          qTerm = gsl_complex_mul_real(qTerm, laser_fe.polarization[q]);
          oCoherences = gsl_complex_add(oCoherences, qTerm);
        }
        oCoherences = gsl_complex_mul_imag(oCoherences,
                                           laser_fe.field/(2*_planck_hbar));
        drho_ff[f][fp] = gsl_complex_add(drho_ff[f][fp], oCoherences);
      } // End if doing coherences
    }   // End fp loop
  }     // End f loop
  // G-Ground states rho_gg:  Equation 34
  // Zeeman coherences are when g != g`
  for (int g = 0; g < numGStates; g++) {
    for (int gp = 0; gp < numGStates; gp++) {
      gsl_complex oCoherences = gsl_complex_rect(0.0, 0.0);
      if ((g == gp) || gs_Zeeman) {
        for (int q = 0; q < 3; q++) {
          gsl_complex qTerm = gsl_complex_rect(0.0, 0.0);
          for (int epp = 0; epp < numEStates; epp++) {
            gsl_complex left, right;  // First and second halves in
            // top line of equation 34
            left = gsl_complex_rect(pow(-1.0, (0))*
                                    dipole_moment_eg[epp][g]*a_eg[epp][g][q],
                                    0.0);
            left = gsl_complex_mul(left,
                                   delta_eg[epp][gp]);
            right = gsl_complex_rect(dipole_moment_eg[epp][gp]*a_eg[epp][gp][q],
                                     0.0);
            right = gsl_complex_mul(right,
                                    gsl_complex_conjugate(delta_eg[epp][g]));
            left = gsl_complex_sub(left, right);
            qTerm = gsl_complex_add(qTerm, left);
          }
          qTerm = gsl_complex_mul_real(qTerm, laser_ge.polarization[q]);
          oCoherences = gsl_complex_add(oCoherences, qTerm);
        }
        oCoherences = gsl_complex_mul_imag(oCoherences,
                                           laser_ge.field/(2*_planck_hbar));
        drho_gg[g][gp] = gsl_complex_add(drho_gg[g][gp], oCoherences);
      } // End if doing coherences
    }   // End gp loop
  }     // End g loop
  }
  for (int e = 0; e < numEStates; e++) {
    for (int ep = 0; ep < numEStates; ep++) {
      drho_ee[e][ep] = gsl_complex_mul_real(drho_ee[e][ep], dt);
      rho_ee[e][ep] = gsl_complex_add(rho_ee[e][ep], drho_ee[e][ep]);
    }
    for (int f = 0; f < numFStates; f++) {
      ddelta_ef[e][f] = gsl_complex_mul_real(ddelta_ef[e][f], dt);
      delta_ef[e][f] = gsl_complex_add(delta_ef[e][f], ddelta_ef[e][f]);
    }
    for (int g = 0; g < numGStates; g++) {
      ddelta_eg[e][g] = gsl_complex_mul_real(ddelta_eg[e][g], dt);
      delta_eg[e][g] = gsl_complex_add(delta_eg[e][g], ddelta_eg[e][g]);
    }
  }
  for (int f = 0; f < numFStates; f++) {
    for (int fp = 0; fp < numFStates; fp++) {
      drho_ff[f][fp] = gsl_complex_mul_real(drho_ff[f][fp], dt);
      rho_ff[f][fp] = gsl_complex_add(rho_ff[f][fp], drho_ff[f][fp]);
    }
  }
  for (int g = 0; g < numGStates; g++) {
    for (int gp = 0; gp < numGStates; gp++) {
      drho_gg[g][gp] = gsl_complex_mul_real(drho_gg[g][gp], dt);
      rho_gg[g][gp] = gsl_complex_add(rho_gg[g][gp], drho_gg[g][gp]);
    }
  }
}

void Density_Matrix::setup_dipole_moments(double gamma) {
  printf("Density_Matrix::setup_dipole_moments(double gamma)\n");
  for (int e = 0; e < numEStates; e++) {
    double ang_freq_ex = nu_E[e] * 2.0 * M_PI;
    for (int g = 0; g < numGStates; g++) {
      printf("e = %d <---> g = %d \t", e, g);
      double ang_freq_gr = nu_G[g] * 2.0 * M_PI;
      dipole_moment_eg[e][g] = set_dipole_moment(gamma, ang_freq_ex,
                                                ang_freq_gr);
    }
    for (int f = 0; f < numFStates; f++) {
      printf("e = %d <---> f = %d\t", e, f);
      double ang_freq_gr = nu_F[f] * 2.0 * M_PI;
      dipole_moment_ef[e][f] = set_dipole_moment(gamma, ang_freq_ex,
                                                 ang_freq_gr);
    }
  }
}

double Density_Matrix::set_dipole_moment(double gamma, double omega_ex,
                                         double omega_gr) {
  // This calculates the dipole moment matrix element given by Tremblay 1990
  // equation 31.  It does not include reference to the polarization (e_q) or
  // reference to the transfert coefficients (a_eg).  Basically this corresponds
  // to the radial part of the matrix element
  double num = 3*M_PI*_epsilon_0*_planck_hbar*pow(_speed_of_light, 3.0) * gamma;
  double den = pow(omega_ex-omega_gr, 3.0);
  double dipole = sqrt(num/den);
  printf("gamma = %5.3G MHz\t delta_omega = %8.6G MHz\t", gamma/_MHz,
         fabs(omega_ex-omega_gr)/_MHz);
  printf("dipole moment = %8.6G e*nm\n", dipole/(_elementary_charge*_nm));
  return dipole;
}

void Density_Matrix::print_rabi_frequencies(FILE * des) {
  // Calculate Rabi frequencies based on Metcalf equation 1.10.  Note that this
  // equation has the angular frequency!  I will also print the period to avoid
  // confusion.
  fprintf(des, "*****Rabi Frequencies*****\n");
  for (int q = 0; q < 3; q++) {
    if (laser_fe.polarization[q] > 0.0 || laser_ge.polarization[q] > 0.0)
      fprintf(des, "***** q = %d *****\n", q-1);
    for (int e = 0; e < numEStates; e++) {
      if (laser_ge.polarization[q] > 0.0) {
        for (int g = 0; g < numGStates; g++) {
          double rabi = dipole_moment_eg[e][g] * a_eg[e][g][q] * laser_ge.field;
          rabi = rabi / _planck_hbar;
          double detune = 2*M_PI*(laser_ge.nu - (nu_E[e] - nu_G[g]));
          rabi = sqrt(pow(rabi, 2.0) + pow(detune, 2.0));
          double period = 2*M_PI / rabi;
          fprintf(des, "e = %d <--> g = %d\t", e, g);
          fprintf(des, "dipole moment = %8.6G e*nm\t",
                  dipole_moment_eg[e][g]/(_elementary_charge*_nm));
          fprintf(des, "detune = %8.6G MHz\t", detune/_MHz);
          fprintf(des, "rabi frequency = %8.6G MHz\t", rabi/_MHz);
          fprintf(des, "period = %8.6G us\n", period/_us);
        }
      }
      if (laser_fe.polarization[q] > 0.0) {
        for (int f = 0; f < numFStates; f++) {
          double rabi = dipole_moment_ef[e][f] * a_ef[e][f][q] * laser_fe.field;
          rabi = rabi / _planck_hbar;
          double detune = 2*M_PI*(laser_fe.nu - (nu_E[e] - nu_F[f]));
          //        rabi = sqrt(pow(rabi, 2.0) + pow(detune, 2.0));
          double period = 2*M_PI / rabi;
          fprintf(des, "e = %d <--> f = %d\t", e, f);
          fprintf(des, "dipole moment = %8.6G e*nm\t",
                  dipole_moment_ef[e][f]/(_elementary_charge*_nm));
          fprintf(des, "detune = %8.6G MHz\t", detune/_MHz);
          fprintf(des, "rabi frequency = %8.6G MHz\t", rabi/_MHz);
          fprintf(des, "period = %8.6G us\n", period/_us);
        }
      }
    }
  }
}

void Density_Matrix::reset_dPop() {
  for (int e = 0; e < numEStates; e++) {
    for (int ep = 0; ep < numEStates; ep++) {
      drho_ee[ep][ep] = gsl_complex_rect(0.0, 0.0);
    }
    for (int f = 0; f < numFStates; f++) {
      ddelta_ef[f][f] = gsl_complex_rect(0.0, 0.0);
    }
    for (int g = 0; g < numGStates; g++) {
      ddelta_eg[g][g] = gsl_complex_rect(0.0, 0.0);
    }
  }
  for (int f = 0; f < numFStates; f++) {
    for (int fp = 0; fp < numFStates; fp++) {
      drho_ff[fp][fp] = gsl_complex_rect(0.0, 0.0);
    }
    for (int g = 0; g < numGStates; g++) {
      ddelta_fg[g][g] = gsl_complex_rect(0.0, 0.0);
    }
  }
  for (int g = 0; g < numGStates; g++) {
    for (int gp = 0; gp < numGStates; gp++) {
      drho_gg[g][g] = gsl_complex_rect(0.0, 0.0);
    }
  }
}
