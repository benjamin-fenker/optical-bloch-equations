// Authors: Benjamin Fenker 2013
// Copyright 2012 Benjamin Fenker
#include <stdio.h>
#include <math.h>
#include "include/rate_equations.h"
#include "include/units.h"

using std::vector;
extern bool op_verbose;

//#pragma omp declare reduction (gsl_complex_add : gsl_complex : omp_out = gsl_complex_add(omp_out, omp_in))

Rate_Equations::Rate_Equations() {
}

Rate_Equations::Rate_Equations(Eigenvector_Helper set_eigen,
                               Laser_data set_laser_fe, Laser_data set_laser_ge,
                               double tilt)
    : OpticalPumping_Method(set_eigen, set_laser_fe, set_laser_ge, tilt),
      transition_rate_eg(eigen.atom.numEStates,
                         vector<vector<double> >(eigen.atom.numGStates,
                                                 vector<double>(3, 0.0))),
      transition_rate_ef(eigen.atom.numEStates,
                         vector<vector<double> >(eigen.atom.numFStates,
                                                 vector<double>(3, 0.0))),
      dPop_g(eigen.atom.numGStates, 0.0), dPop_f(eigen.atom.numFStates, 0.0),
      dPop_e(eigen.atom.numEStates, 0.0) {
  //   printf("Rate_Equations::Rate_Equations(...)\n\n");
  if (op_verbose) {
    printf("Stokes vector: <%8.6G, %4.2G, %4.2G, %8.6G\n", laser_fe.stokes[0],
           laser_fe.stokes[1], laser_fe.stokes[2], laser_fe.stokes[3]);
  }
  
  setup_transition_rates(eigen.atom.linewidth);
  totalTerms = eigen.atom.numFStates+eigen.atom.numGStates;
  totalTerms += eigen.atom.numEStates;
  es_Zeeman = false;
  gs_Zeeman = false;
  es_hyperfine = false;
  gs_hyperfine = false;
}

Rate_Equations::~Rate_Equations() {
  // delete population;
}

void Rate_Equations::calculate_derivs(DM_container *status) {
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
          double pop_diff = (GSL_REAL(status->gg[g][g])-
                             GSL_REAL(status->ee[e][e]));
          if (debug) printf("|g=%d> <--> |e%d\t", g, e);
          if (debug) printf(" popDiff = %8.6G\t tran_rate = %8.6G ns^-1\t",
                            pop_diff, transition_rate_eg[e][g][q]/(1/_ns));
          if (debug) printf("a^2 = %5.3G\t", pow(a_eg[e][g][q], 2.0));
          delta_laser = pop_diff * transition_rate_eg[e][g][q];
          delta_laser *= pow(a_eg[e][g][q], 2.0);
          if (debug) printf("delta_laser = %8.6G ns^-1\n", delta_laser/(1/_ns));
          dm_derivs->ee[e][e] = gsl_complex_add_real(dm_derivs->ee[e][e],
                                                     delta_laser);
          dm_derivs->gg[g][g] = gsl_complex_sub_real(dm_derivs->gg[g][g],
                                                     delta_laser);
        // Spontaneous emission
        delta_spon = GSL_REAL(status->ee[e][e]) *
            pow(a_eg[e][g][q], 2.0) / tau;
        dm_derivs->ee[e][e] = gsl_complex_sub_real(dm_derivs->ee[e][e],
                                                   delta_spon);
        dm_derivs->gg[g][g] = gsl_complex_add_real(dm_derivs->gg[g][g],
                                                   delta_spon);
        }  // End q-loop
      }  // End g-loop
      // f <---> e
      for (int f =0; f < numFStates; f++) {
        for (int q = 0; q < 3; q++) {
          double pop_diff = (GSL_REAL(status->ff[f][f])-
                             GSL_REAL(status->ee[e][e]));
          if (debug) printf("|f=%d> <--> |e=%d\t ", f, e);
          if (debug) printf("popDiff = %8.6G\t tran_rate = %8.6G ns^-1\t",
                            pop_diff, transition_rate_ef[e][f][q]/(1/_ns));
          if (debug) printf("a^2 = %5.3G\t", pow(a_ef[e][f][q], 2.0));
          delta_laser = pop_diff * transition_rate_ef[e][f][q];
          delta_laser *=  pow(a_ef[e][f][q], 2.0);
          if (debug) printf("delta_laser = %8.6G ns^-1\n", delta_laser/(1/_ns));
          dm_derivs->ee[e][e] = gsl_complex_add_real(dm_derivs->ee[e][e],
                                                     delta_laser);
          dm_derivs->ff[f][f] = gsl_complex_sub_real(dm_derivs->ff[f][f],
                                                     delta_laser);
          // Spontaneous emission
          delta_spon = GSL_REAL(status->ee[e][e]) *
              pow(a_ef[e][f][q], 2.0) / tau;
          dm_derivs->ee[e][e] = gsl_complex_sub_real(dm_derivs->ee[e][e],
                                                     delta_spon);
          dm_derivs->ff[f][f] = gsl_complex_add_real(dm_derivs->ff[f][f],
                                                     delta_spon);
        }  // End q-loop
      }  // End f-loop
    }    // End e-loop
    //  apply_transverse_field(status);
}

void Rate_Equations::setup_transition_rates(double linewidth) {
  bool debug = true;
  if (debug) {
    printf("Laser_g %8.6G MHz\t Laser_f %8.6G MHz\n", laser_ge.nu/_MHz,
           laser_fe.nu/_MHz);
  }
  for (int q = 0; q < 3; q++) {
    for (int e = 0; e < numEStates; e++) {
      for (int g = 0; g < numGStates; g++) {
        double transition_freq = nu_E[e] - nu_G[g];
        if (debug) {
          printf("|g=%d>-->|e=%d (%d)\t nu_eg = %14.10G MHz\t", g, e, q,
                 transition_freq/_MHz);
          printf("nu_L (ge) = %14.10G MHz\n", laser_ge.nu/_MHz);
        }
        // The tuned laser
        transition_rate_eg[e][g][q] = set_transition_rate(
            laser_ge.intensity[q], laser_ge.saturation_intensity, linewidth,
            laser_ge.linewidth, transition_freq, laser_ge.nu);
        // The other (off resonance) laser
        transition_rate_eg[e][g][q] += set_transition_rate(
            laser_fe.intensity[q], laser_fe.saturation_intensity, linewidth,
            laser_fe.linewidth, transition_freq, laser_fe.nu);
      }  // End g-loop
      for (int f = 0; f < numFStates; f++) {
        double transition_freq = nu_E[e] - nu_F[f];
        if (debug) {
          printf("|f=%d>-->|e=%d (%d)\t nu_eg = %14.10G MHz\t", f, e, q,
                 transition_freq);
          printf("nu_L (fe) = %14.10G MHz\n", laser_fe.nu/_MHz);
        }
        transition_rate_ef[e][f][q] = set_transition_rate(
            laser_fe.intensity[q], laser_fe.saturation_intensity, linewidth,
            laser_fe.linewidth, transition_freq, laser_fe.nu);
        transition_rate_ef[e][f][q] += set_transition_rate(
            laser_ge.intensity[q], laser_ge.saturation_intensity, linewidth,
            laser_ge.linewidth, transition_freq, laser_ge.nu);
      }  // End f-loop
    }    // End e-loop
  }      // End q-loop
}        // End setup

double Rate_Equations::set_transition_rate(double laser_power,
                                           double sat_intensity,
                                           double atom_lw,
                                           double laser_lw, double atom_freq,
                                           double laser_freq) {
  bool debug = true;
  // This uses Nafcha equation 15 and Metcalf equation 2.24c to define the
  // saturation intensity
  // Note that Nafcha's linewidths of FWHM/2.  My linewidths are defined as the
  // FWHM of the laser
  double rate = laser_power / (4.0 * M_PI * pow(tau, 2.0) * sat_intensity);
  double lorentzian = 4.0*pow(laser_freq - atom_freq, 2.0) +
    pow(laser_lw + atom_lw, 2.0);
  //  lorentzian = pow(laser_lw + atom_lw, 2.0);
  // Uncommenting the line above turns off the detuning factor
  lorentzian = (laser_lw + atom_lw)/lorentzian;
  rate *= lorentzian;                   // The right way

  // DELETE ME 
  //  if (rate < 0.5 * _MHz) rate = 0.0;
  // DELETE ME

  if (debug && rate > 0.0) {
    printf("\tDetune = %8.6G MHz\tLorentzian = %8.6G ns\t rate = %10.8G MHz\n",
           ((laser_freq-atom_freq))/_MHz, lorentzian/_ns, rate/_MHz);
  }

  return rate;
}

void Rate_Equations::switch_off_laser(int las) {
  OpticalPumping_Method::switch_off_laser(las);
  setup_transition_rates(eigen.atom.linewidth);
}

void Rate_Equations::change_magnetic_field(double newField) {
  OpticalPumping_Method::change_magnetic_field(newField);
  setup_transition_rates(eigen.atom.linewidth);
}
