// Copyright 2012 Benjamin Fenker
#include <stdio.h>
#include <math.h>
#include "include/rate_equations.h"
#include "include/units.h"

using std::vector;

Rate_Equations::Rate_Equations(atom_data atom, magnetic_field_data field,
                               Laser_data set_laser_fe, Laser_data set_laser_ge)
  : OpticalPumping_Method(atom, field, set_laser_fe, set_laser_ge),
    transition_rate_eg(atom.numEStates, vector<double>(atom.numGStates, 0.0)),
    transition_rate_ef(atom.numEStates, vector<double>(atom.numFStates, 0.0)),
    dPop_g(atom.numGStates, 0.0), dPop_f(atom.numFStates, 0.0),
    dPop_e(atom.numEStates, 0.0) {
  printf("Rate_Equations::Rate_Equations(...)\n");
  print_couplings(stdout);
  setup_transition_rates(atom.linewidth);

  totalTerms = atom.numFStates+atom.numGStates+atom.numEStates;


  data.gStart = 0;
  data.gEnd = data.numGStates;   // Means I'm going to use g < gEnd rather than
  // g <= gEnd
  data.fStart = data.gEnd;
  data.fEnd = data.fStart + data.numFStates;
  data.eStart = data.fEnd;
  data.eEnd = data.eStart + data.numEStates;
  data.transition_rate_eg = transition_rate_eg;
  data.transition_rate_ef = transition_rate_ef;

  population = new double[atom.numFStates+atom.numGStates+atom.numEStates];
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

int Rate_Equations::update_population_gsl(double t, const double y[],
                                                 double f[], void *params) {
  op_data_for_gsl data = *(reinterpret_cast<op_data_for_gsl *>(params));
  //  printf("t = %5.3G   tau = %5.3G ns\t", t, data.atom.tau/_ns);

  // Markers to keep track of where I am in the array!
  //  printf("g: %d --> %d\tf: %d --> %d\t e: %d-->%d\n", gStart, gEnd,
  //      fStart, fEnd, eStart, eEnd);
  for (int i = data.gStart; i < data.eEnd; i++) f[i] = 0.0;
  // This follows exactly what is done in Nafcha 1995 and DM's thesis.
  // I use the expressions from Nafcha as they are more clear where frequency
  // and angular frequency are concerned
  bool debug = false;
  double delta_laser, delta_spon;
  // Laser interaction:  stimulated emission and absoprtion
  for (int q = 0; q < 3; q++) {
    for (int e = data.eStart; e < data.eEnd; e++) {
      int e_eff = e - data.eStart;
      // g <---> e
      for (int g = data.gStart; g < data.gEnd; g++) {
        int g_eff = g - data.gStart;
        if (data.laser_ge.polarization[q] > 0.0) {
          double pop_diff = y[g] - y[e];
          if (debug) printf("|g=%d> <--> |e%d\t", g_eff, e_eff);
          if (debug) printf(" popdiff = %8.6G\t tran_rate = %8.6G ns^-1\t",
                            pop_diff,
                            data.transition_rate_eg[e_eff][g_eff]/(1/_ns));
          if (debug) printf("a^2 = %5.3G\t polarization = %5.3G\t",
                            pow(data.a_eg[e_eff][g_eff][q], 2.0),
                            data.laser_ge.polarization[q]);
          delta_laser = pop_diff * data.transition_rate_eg[e_eff][g_eff];
          delta_laser *= pow(data.a_eg[e_eff][g_eff][q], 2.0);
          delta_laser *= data.laser_ge.polarization[q];
          if (debug) printf("delta_laser = %8.6G ns^-1\t", delta_laser/(1/_ns));
          f[e] += delta_laser;
          f[g] -= delta_laser;
        }
        // Spontaneous emission
        delta_spon = y[e] * pow(data.a_eg[e_eff][g_eff][q], 2.0);
        delta_spon /= data.atom.tau;
        f[e] -= delta_spon;
        f[g] += delta_spon;
      }
      // f <---> e
      for (int fi = data.fStart; fi < data.fEnd; fi++) {
        int f_eff = fi - data.fStart;
        if (data.laser_fe.polarization[q] > 0.0) {
          double pop_diff = y[fi] - y[e];
          if (debug) printf("|f=%d> <--> |e=%d\t ", fi, e);
          if (debug) printf("popDiff = %8.6G\t tran_rate = %8.6G ns^-1\t",
                            pop_diff,
                            data.transition_rate_ef[e_eff][f_eff]/(1/_ns));
          if (debug) printf("a^2 = %5.3G\t polarization = %5.3G\t",
                            pow(data.a_ef[e_eff][f_eff][q], 2.0),
                            data.laser_fe.polarization[q]);
          delta_laser = pop_diff * data.transition_rate_ef[e_eff][f_eff];
          delta_laser *=  pow(data.a_ef[e_eff][f_eff][q], 2.0);
          delta_laser *=  data.laser_fe.polarization[q];
          if (debug) printf("delta_laser = %8.6G ns^-1\t", delta_laser/(1/_ns));
          f[e] += delta_laser;
          f[fi] -= delta_laser;
        }
        // Spontaneous emission
        delta_spon = y[e] * pow(data.a_ef[e_eff][f_eff][q], 2.0);
        delta_spon /= data.atom.tau;
        f[e] -= delta_spon;
        f[fi] += delta_spon;
      }
    }
  }
  return GSL_SUCCESS;
}


void Rate_Equations::update_population(double dt) {
  // This follows exactly what is done in Nafcha 1995 and DM's thesis.
  // I use the expressions from Nafcha as they are more clear where frequency
  // and angular frequency are concerned
  bool debug = false;
  reset_dPop();
  double delta_laser, delta_spon;
  // Laser interaction:  stimulated emission and absoprtion
  for (int q = 0; q < 3; q++) {
    for (int e = 0; e < numEStates; e++) {
      // g <---> e
      for (int g =0; g < numGStates; g++) {
        if (laser_ge.polarization[q] > 0.0) {
          double pop_diff = (GSL_REAL(rho_gg[g][g])-GSL_REAL(rho_ee[e][e]));
          if (debug) printf("|g=%d> <--> |e%d\t", g, e);
          if (debug) printf(" popDiff = %8.6G\t tran_rate = %8.6G ns^-1\t",
                            pop_diff, transition_rate_eg[e][g]/(1/_ns));
          if (debug) printf("a^2 = %5.3G\t polarization = %5.3G\t",
                            pow(a_eg[e][g][q], 2.0), laser_ge.polarization[q]);
          delta_laser = pop_diff * transition_rate_eg[e][g];
          delta_laser *= pow(a_eg[e][g][q], 2.0) * laser_ge.polarization[q];
          if (debug) printf("delta_laser = %8.6G ns^-1\t", delta_laser/(1/_ns));
          dPop_e[e] += delta_laser;
          dPop_g[g] -= delta_laser;
          if (debug) printf("dPop_e[%d] = %12.10G\t dPop_g[%d] = %12.10G\n",
                            e, dt*dPop_e[e], g, dt*dPop_g[g]);
        }
        // Spontaneous emission
        delta_spon = GSL_REAL(rho_ee[e][e]) * pow(a_eg[e][g][q], 2.0) / tau;
        dPop_e[e] -= delta_spon;
        dPop_g[g] += delta_spon;
      }
      // f <---> e
      for (int f =0; f < numFStates; f++) {
        if (laser_fe.polarization[q] > 0.0) {
          double pop_diff = (GSL_REAL(rho_ff[f][f])-GSL_REAL(rho_ee[e][e]));
          if (debug) printf("|f=%d> <--> |e=%d\t ", f, e);
          if (debug) printf("popDiff = %8.6G\t tran_rate = %8.6G ns^-1\t",
                            pop_diff, transition_rate_ef[e][f]/(1/_ns));
          if (debug) printf("a^2 = %5.3G\t polarization = %5.3G\t",
                            pow(a_ef[e][f][q], 2.0),
                            laser_fe.polarization[q]);
          delta_laser = pop_diff * transition_rate_ef[e][f];
          delta_laser *=  pow(a_ef[e][f][q], 2.0) * laser_fe.polarization[q];
          if (debug) printf("delta_laser = %8.6G ns^-1\t", delta_laser/(1/_ns));
          dPop_e[e] += delta_laser;
          dPop_f[f] -= delta_laser;
          if (debug) printf("dPop_e[%d] = %12.10G\t dPop_f[%d] = %12.10G\n",
                            e, dt*dPop_e[e], f, dt*dPop_f[f]);
        }
        // Spontaneous emission
        delta_spon = GSL_REAL(rho_ee[e][e]) * pow(a_ef[e][f][q], 2.0) / tau;
        dPop_e[e] -= delta_spon;
        dPop_f[f] += delta_spon;
      }
    }
  }

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
  printf("G = %d F = %d E = %d\n", numGStates, numFStates, numEStates);
  printf("Laser_g %8.6G MHz\t Laser_f %8.6G MHz\n", laser_ge.nu/_MHz,
         laser_fe.nu/_MHz);
  for (int e = 0; e < numEStates; e++) {
    for (int g = 0; g < numGStates; g++) {
      double transition_freq = nu_E[e] - nu_G[g];
      printf("|g=%d>-->|e=%d\t nu_eg = %14.10G MHz\tnu_L = %14.10G MHz\t", g, e,
             transition_freq/_MHz, laser_ge.nu/_MHz);
      transition_rate_eg[e][g] = set_transition_rate(
                                       laser_ge.power,
                                        laser_ge.saturation_intensity,
                                        linewidth, laser_ge.linewidth,
                                        transition_freq, laser_ge.nu);
    }
    for (int f = 0; f < numFStates; f++) {
      double transition_freq = nu_E[e] - nu_F[f];
      printf("|f=%d>-->|e=%d\t nu_eg = %14.10G MHz\tnu_L = %14.10G MHz\t", f, e,
             transition_freq/_MHz, laser_fe.nu/_MHz);
      transition_rate_ef[e][f] = set_transition_rate(
                                     laser_fe.power,
                                     laser_fe.saturation_intensity,
                                     linewidth, laser_fe.linewidth,
                                     transition_freq, laser_fe.nu);
    }
  }
}

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
  double lorentzian = pow(laser_freq- atom_freq, 2.0) +
    4.0 * pow(laser_lw + atom_lw, 2.0);
  // lorentzian = pow(laser_lw + atom_lw, 2.0);
  // Uncommenting the line above turns off the detuning factor
  lorentzian = (laser_lw + atom_lw)/lorentzian;
  rate *= lorentzian;
  printf("Detune = %8.6G MHz\tLorentzian = %8.6G ns\t rate = %8.6G MHz\n",
         (fabs(laser_freq-atom_freq))/_MHz, lorentzian/_ns, rate/_MHz);
  return rate;
}
