// Authors: Benjamin Fenker 2013
// Copyright 2012 Benjamin Fenker
#include <stdio.h>
#include <math.h>
#include "include/rate_equations.h"
#include "include/units.h"

using std::vector;
extern bool op_verbose;

Rate_Equations::Rate_Equations(Eigenvector_Helper set_eigen,
                               Laser_data set_laser_fe, Laser_data set_laser_ge)
  : OpticalPumping_Method(set_eigen, set_laser_fe, set_laser_ge),
    transition_rate_eg(eigen.atom.numEStates,
                       vector<vector<double> >(eigen.atom.numGStates,
                                               vector<double>(3, 0.0))),
    transition_rate_ef(eigen.atom.numEStates,
                       vector<vector<double> >(eigen.atom.numFStates,
                                               vector<double>(3, 0.0))),
    dPop_g(eigen.atom.numGStates, 0.0), dPop_f(eigen.atom.numFStates, 0.0),
    dPop_e(eigen.atom.numEStates, 0.0) {
  printf("Rate_Equations::Rate_Equations(...)\n\n");
  if (op_verbose) {
    printf("Stokes vector: <%8.6G, %4.2G, %4.2G, %8.6G\n", laser_fe.stokes[0],
           laser_fe.stokes[1], laser_fe.stokes[2], laser_fe.stokes[3]);
  }

  setup_transition_rates(eigen.atom.linewidth);
  totalTerms = eigen.atom.numFStates+eigen.atom.numGStates;
  totalTerms += eigen.atom.numEStates;

  /* 
  data.gStart = 0;
  data.gEnd = data.numGStates;   // Means I'm going to use g < gEnd rather than
  // g <= gEnd
  data.fStart = data.gEnd;
  data.fEnd = data.fStart + data.numFStates;
  data.eStart = data.fEnd;
  data.eEnd = data.eStart + data.numEStates;
  data.transition_rate_eg = transition_rate_eg;
  data.transition_rate_ef = transition_rate_ef;
  */

  population = new double[eigen.atom.numFStates+eigen.atom.numGStates+
                          eigen.atom.numEStates];
  int i, j;
  for (i = 0; i < (numFStates+numGStates); i++) {
    population[i] = 1.0 / (numFStates + numGStates);
  }
  for (j = i; j < totalTerms; j++) {
    population[j] = 0.0;
  }
}

Rate_Equations::~Rate_Equations() {
  delete population;
}

void Rate_Equations::update_population(double dt) {
  // This follows exactly what is done in Nafcha 1995 and DM's thesis.
  // I use the expressions from Nafcha as they are more clear where frequency
  // and angular frequency are concerned
  bool debug = false;
  reset_dPop();
  //  printf("RATE EQUATION!!!\n");
  double delta_laser, delta_spon;
  // Laser interaction:  stimulated emission and absoprtion
    for (int e = 0; e < numEStates; e++) {
      // g <---> e
      for (int g = 0; g < numGStates; g++) {
        for (int q = 0; q < 3; q++) {
          double pop_diff = (GSL_REAL(rho_gg[g][g])-GSL_REAL(rho_ee[e][e]));
          if (debug) printf("|g=%d> <--> |e%d\t", g, e);
          if (debug) printf(" popDiff = %8.6G\t tran_rate = %8.6G ns^-1\t",
                            pop_diff, transition_rate_eg[e][g][q]/(1/_ns));
          if (debug) printf("a^2 = %5.3G\t", pow(a_eg[e][g][q], 2.0));
          delta_laser = pop_diff * transition_rate_eg[e][g][q];
          delta_laser *= pow(a_eg[e][g][q], 2.0);
          if (debug) printf("delta_laser = %8.6G ns^-1\t", delta_laser/(1/_ns));
          dPop_e[e] += delta_laser;
          dPop_g[g] -= delta_laser;
          if (debug) printf("dPop_e[%d] = %12.10G\t dPop_g[%d] = %12.10G\n",
                            e, dt*dPop_e[e], g, dt*dPop_g[g]);
        // Spontaneous emission
        delta_spon = GSL_REAL(rho_ee[e][e]) * pow(a_eg[e][g][q], 2.0) / tau;
        dPop_e[e] -= delta_spon;
        dPop_g[g] += delta_spon;
        }  // End q-loop
      }  // End g-loop
      // f <---> e
      for (int f =0; f < numFStates; f++) {
        for (int q = 0; q < 3; q++) {
          double pop_diff = (GSL_REAL(rho_ff[f][f])-GSL_REAL(rho_ee[e][e]));
          if (debug) printf("|f=%d> <--> |e=%d\t ", f, e);
          if (debug) printf("popDiff = %8.6G\t tran_rate = %8.6G ns^-1\t",
                            pop_diff, transition_rate_ef[e][f][q]/(1/_ns));
          if (debug) printf("a^2 = %5.3G\t", pow(a_ef[e][f][q], 2.0));
          delta_laser = pop_diff * transition_rate_ef[e][f][q];
          delta_laser *=  pow(a_ef[e][f][q], 2.0);
          if (debug) printf("delta_laser = %8.6G ns^-1\t", delta_laser/(1/_ns));
          dPop_e[e] += delta_laser;
          dPop_f[f] -= delta_laser;
          if (debug) printf("dPop_e[%d] = %12.10G\t dPop_f[%d] = %12.10G\n",
                            e, dt*dPop_e[e], f, dt*dPop_f[f]);
          // Spontaneous emission
          delta_spon = GSL_REAL(rho_ee[e][e]) * pow(a_ef[e][f][q], 2.0) / tau;
          dPop_e[e] -= delta_spon;
          dPop_f[f] += delta_spon;
        }  // End q-loop
      }  // End f-loop
    }    // End e-loop

  for (int g = 0; g < numGStates; g++) {
    GSL_SET_REAL(&rho_gg[g][g], GSL_REAL(rho_gg[g][g]) + dt*dPop_g[g]);
  }
  for (int e = 0; e < numEStates; e++) {
    GSL_SET_REAL(&rho_ee[e][e], GSL_REAL(rho_ee[e][e]) + dt*dPop_e[e]);
  }
  for (int f = 0; f < numFStates; f++) {
    GSL_SET_REAL(&rho_ff[f][f], GSL_REAL(rho_ff[f][f]) + dt*dPop_f[f]);
  }
}

void Rate_Equations::reset_dPop() {
  for (int g = 0; g < numGStates; g++) dPop_g[g] = 0.0;
  for (int f = 0; f < numFStates; f++) dPop_f[f] = 0.0;
  for (int e = 0; e < numEStates; e++) dPop_e[e] = 0.0;
}

void Rate_Equations::setup_transition_rates(double linewidth) {
  if (op_verbose) {
    printf("Laser_g %8.6G MHz\t Laser_f %8.6G MHz\n", laser_ge.nu/_MHz,
           laser_fe.nu/_MHz);
  }
  for (int q = 0; q < 3; q++) {
    for (int e = 0; e < numEStates; e++) {
      for (int g = 0; g < numGStates; g++) {
        double transition_freq = nu_E[e] - nu_G[g];
        if (op_verbose) {
          printf("|g=%d>-->|e=%d (%d)\t nu_eg = %14.10G MHz\t", g, e, q,
                 transition_freq/_MHz);
          printf("nu_L = %14.10G MHz\t", laser_ge.nu/_MHz);
        }
        transition_rate_eg[e][g][q] = set_transition_rate(
                                          laser_ge.intensity[q],
                                          laser_ge.saturation_intensity,
                                          linewidth, laser_ge.linewidth,
                                          transition_freq, laser_ge.nu);
      }  // End g-loop
      for (int f = 0; f < numFStates; f++) {
        double transition_freq = nu_E[e] - nu_F[f];
        if (op_verbose) {
          printf("|f=%d>-->|e=%d (%d)\t nu_eg = %14.10G MHz\t", f, e, q,
                 transition_freq);
          printf("nu_L = %14.10G MHz\t", laser_fe.nu/_MHz);
        }
        transition_rate_ef[e][f][q] = set_transition_rate(
                                          laser_fe.intensity[q],
                                          laser_fe.saturation_intensity,
                                          linewidth, laser_fe.linewidth,
                                          transition_freq, laser_fe.nu);
      }  // End f-loop
    }    // End e-loop
  }      // End q-loop
}        // End setup

double Rate_Equations::set_transition_rate(double laser_power,
                                           double sat_intensity,
                                           double atom_lw,
                                           double laser_lw, double atom_freq,
                                           double laser_freq) {
  // This uses Nafcha equation 15 and Metcalf equation 2.24c to define the
  // saturation intensity
  // Note that Nafcha's linewidths of FWHM/2.  My linewidths are defined as the
  // FWHM of the laser
  double rate = laser_power / (4.0 * M_PI * pow(tau, 2.0) * sat_intensity);
  double lorentzian = 4.0*pow(laser_freq - atom_freq, 2.0) +
    pow(laser_lw + atom_lw, 2.0);
  // lorentzian = pow(laser_lw + atom_lw, 2.0);
  // Uncommenting the line above turns off the detuning factor
  lorentzian = (laser_lw + atom_lw)/lorentzian;
  rate *= lorentzian;
  if (op_verbose) {
    printf("Detune = %8.6G MHz\tLorentzian = %8.6G ns\t rate = %8.6G MHz\n",
           (fabs(laser_freq-atom_freq))/_MHz, lorentzian/_ns, rate/_MHz);
  }
  return rate;
}
