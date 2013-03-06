// Authors: Benjamin Fenker 2013

// Copyright Benjamin Fenker 2013
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include "include/density_matrix.h"
#include "include/units.h"

using std::vector;
using std::min;
extern bool op_verbose;

Density_Matrix::Density_Matrix(Eigenvector_Helper set_eigen,
                               Laser_data set_laser_fe, Laser_data set_laser_ge,
                               coherence_flags flags)
  : OpticalPumping_Method(set_eigen, set_laser_fe, set_laser_ge),
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
  // printf("DensityMatrix::Density_Matrix(...)\n\n");
  es_Zeeman = flags.zCoherences;
  gs_Zeeman = flags.zCoherences;
  es_hyperfine = flags.hfCoherences_ex;
  gs_hyperfine = flags.hfCoherences_gr;
  setup_dipole_moments(1.0/eigen.atom.tau);
  if (op_verbose) print_rabi_frequencies(stdout);
}

/*
int Density_Matrix::update_population_gsl(double t, const double y[],
                                           double f[], void *params) {
  double mu = *(double *)params;
  printf("t = %4.2G\n", t);
  mu = 10.0;
  f[0] = y[1];
  f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
  return GSL_SUCCESS;

}
*/

void Density_Matrix::update_population(double dt) {
  reset_dPop();
  //  printf("DENSITY MATRIX!!!\n");
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
  for (int e = 0; e < numEStates; e++) {
    for (int ep = 0; ep < numEStates; ep++) {
      if ((e == ep) ||
          (es_hyperfine && MFe2_Vector[e] == MFe2_Vector[ep]) ||
          (es_Zeeman && Fe2_Vector[e] == Fe2_Vector[ep]) ||
          (es_hyperfine && es_Zeeman)) {
        { // G-laser term
          gsl_complex q_sum = gsl_complex_rect(0.0, 0.0);
          for (int q = 0; q < 3; q +=2) {  // Only q = 0, 2 will have lasers
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
              // The 2-q is added due to de Clerq equation I.2
              right = gsl_complex_mul(right,
                                      delta_eg[e][gpp]);
              left = gsl_complex_sub(left, right);
              g_sum = gsl_complex_add(g_sum, left);
            }
            g_sum = gsl_complex_mul_real(g_sum, laser_ge.field[q]);
            q_sum = gsl_complex_add(q_sum, g_sum);
          }
          q_sum = gsl_complex_mul_imag(q_sum, 0.5/ _planck_hbar);
          drho_ee[e][ep] = gsl_complex_add(drho_ee[e][ep], q_sum);
        }  // End G-laser term
        {  // F-laser term
          gsl_complex q_sum = gsl_complex_rect(0.0, 0.0);
          for (int q = 0; q < 3; q += 2) {  // Only q = 0, 2 will have lasers
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
              // The 2-q is added due to de Clerq equation I.2
              right = gsl_complex_mul(right,
                                      delta_ef[e][fpp]);
              left = gsl_complex_sub(left, right);
              f_sum = gsl_complex_add(f_sum, left);
            }
            f_sum = gsl_complex_mul_real(f_sum, laser_fe.field[q]);
            q_sum = gsl_complex_add(q_sum, f_sum);
          }
          q_sum = gsl_complex_mul_imag(q_sum, 0.5/ _planck_hbar);
          drho_ee[e][ep] = gsl_complex_add(drho_ee[e][ep], q_sum);
        }  // End F-Laser term (no loop or if, just an arbitrary block
        {  // Spontaneous decay term
          double angFrequency = 2*M_PI*(nu_E[e] - nu_E[ep]);
          gsl_complex sp_decay = gsl_complex_rect((1.0/data.atom.tau),
                                                  angFrequency);
          sp_decay = gsl_complex_mul(sp_decay, rho_ee[e][ep]);
          drho_ee[e][ep] = gsl_complex_sub(drho_ee[e][ep], sp_decay);
          // printf("Without B_x, drho_ee[%d][%d] = %g + %g i\t", e, ep,
          //        GSL_REAL(drho_ee[e][ep]), GSL_IMAG(drho_ee[e][ep]));
        }  // End spontaneous decay term

        // Get ready for transverse magnetic field!
        gsl_complex Bx = gsl_complex_rect(0.0, 0.0);
        gsl_complex left = gsl_complex_rect(0.0, 0.0);
        gsl_complex rigt = gsl_complex_rect(0.0, 0.0);
        // printf("For e = %d, comparing %d to %d\n", e, MFe2_Vector[e],
        //        MFe2_Vector[e+1]);
        // Is there a state with Mf` = Mf+1 and F` = F?
        if (e+1 < numEStates) {         // Keep things in bounds
          if (MFe2_Vector[e]+2 == MFe2_Vector[e+1]) {  // Sanity check
            gsl_complex tmp = rho_ee[e+1][ep];
            // printf("Bx = %g + %g i + ", GSL_REAL(tmp), GSL_IMAG(tmp));
            gsl_complex_mul_real(tmp, cPlus_E[e]);
            left = gsl_complex_add(left, tmp);
          }
        }
        // Is there a state with Mf` = Mf-1 and F` = F?
        if (e > 0) {                    // Keep things in bounds
          if (MFe2_Vector[e]-2 == MFe2_Vector[e-1]) {  // Sanity check
            gsl_complex tmp = rho_ee[e-1][ep];
            // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
            gsl_complex_mul_real(tmp, cMins_E[e]);
            left = gsl_complex_add(left, tmp);
          }
        }
        left = gsl_complex_mul_real(left, gFactor_E[e]);
        // Is there a state with Mf` = Mf+1 and F` = F?
        if (ep+1 < numEStates) {         // Keep things in bounds
          if (MFe2_Vector[ep]+2 == MFe2_Vector[ep+1]) {  // Sanity check
            gsl_complex tmp = rho_ee[e][ep+1];
            // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
            gsl_complex_mul_real(tmp, cPlus_E[ep]);
            rigt = gsl_complex_sub(rigt, tmp);
          }
        }
        // Is there a state with Mf` = Mf-1 and F` = F?
        if (ep > 0) {         // Keep things in bounds
          if (MFe2_Vector[ep]-2 == MFe2_Vector[ep-1]) {
            gsl_complex tmp = rho_ee[e][ep-1];
            // printf("%g + %g i\n", GSL_REAL(tmp), GSL_IMAG(tmp));
            gsl_complex_mul_real(tmp, cMins_E[ep]);
            rigt = gsl_complex_sub(rigt, tmp);
          }
        }
        rigt = gsl_complex_mul_real(rigt, gFactor_E[ep]);
        Bx = gsl_complex_add(left, rigt);
        Bx = gsl_complex_mul_imag(Bx,
                                  -_bohr_magneton*eigen.field.B_x
                                  /2.0/_planck_hbar);
        // printf("Bx contribution is %g + %g i\n", GSL_REAL(Bx), GSL_IMAG(Bx));
        drho_ee[e][ep] = gsl_complex_add(drho_ee[e][ep], Bx);
      }   // End if coherences
    }      // End ep
  }        // End e

  // G-Ground states rho_gg:  Equation 34
  // Zeeman coherences are when g != g`
  for (int g = 0; g < numGStates; g++) {
    for (int gp = 0; gp < numGStates; gp++) {
      gsl_complex oCoherences = gsl_complex_rect(0.0, 0.0);
      if ((g == gp) || (gs_Zeeman)) {
        for (int q = 0; q < 3; q +=2) {  // Only q = 0, 2 will have lasers
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
            // The 2-q is added due to de Clerq equation I.2
            right = gsl_complex_mul(right,
                                    gsl_complex_conjugate(delta_eg[epp][g]));
            left = gsl_complex_sub(left, right);
            qTerm = gsl_complex_add(qTerm, left);
          }  // End epp loop
          qTerm = gsl_complex_mul_real(qTerm, laser_ge.field[q]);
          oCoherences = gsl_complex_add(oCoherences, qTerm);
        }  // End q-loop
        oCoherences = gsl_complex_mul_imag(oCoherences,
                                           1.0/(2*_planck_hbar));
        drho_gg[g][gp] = gsl_complex_add(drho_gg[g][gp], oCoherences);
        { // Spontaneous decay term
          gsl_complex spon_decay = gsl_complex_rect(0.0, 0.0);
          for (int e = 0; e < numEStates; e++) {
            for (int ep = 0; ep < numEStates; ep++) {
              if (MFe2_Vector[e] - MFe2_Vector[ep] ==
                  MFg2_Vector[g] - MFg2_Vector[gp]) {
                int q = MFe2_Vector[e] - MFg2_Vector[g];
                if (q%2 != 0) {
                  printf("G - spontaneous decay: q2 not even\n");
                  exit(1);
                }
                q = (q/2) + 1;    // Converts q from twice the difference in M_f
                // to an appropriate index
                double coupling = a_eg[e][g][q] * a_eg[ep][gp][q];
                gsl_complex left = gsl_complex_mul_real(rho_ee[e][ep],
                                                        coupling);
                spon_decay = gsl_complex_add(spon_decay, left);
              }  // End if M_f sublevels satisfy m_e - m_e` = m_f - m_f`
            }    // End ep sum
          }      // End e sum
          spon_decay = gsl_complex_mul_real(spon_decay, ((1.0/data.atom.tau)));
          double angFrequency = 2*M_PI*(nu_G[gp] - nu_G[g]);
          gsl_complex right = gsl_complex_mul_imag(rho_gg[g][gp],
                                                   angFrequency);
          spon_decay = gsl_complex_sub(spon_decay, right);
          drho_gg[g][gp] = gsl_complex_add(drho_gg[g][gp], spon_decay);
        }  // End spontaneous decay term
        // Get ready for transverse magnetic field!
        gsl_complex Bx = gsl_complex_rect(0.0, 0.0);
        gsl_complex left = gsl_complex_rect(0.0, 0.0);
        gsl_complex rigt = gsl_complex_rect(0.0, 0.0);
        // printf("For e = %d, comparing %d to %d\n", e, MFe2_Vector[e],
        //        MFe2_Vector[e+1]);
        // Is there a state with Mf` = Mf+1 and F` = F?
        if (g+1 < numGStates) {         // Keep things in bounds
          if (MFg2_Vector[g]+2 == MFg2_Vector[g+1]) {  // Sanity check
            gsl_complex tmp = rho_gg[g+1][gp];
            // printf("Bx = %g + %g i + ", GSL_REAL(tmp), GSL_IMAG(tmp));
            gsl_complex_mul_real(tmp, cPlus_G[g]);
            left = gsl_complex_add(left, tmp);
          }
        }
        // Is there a state with Mf` = Mf-1 and F` = F?
        if (g > 0) {         // Keep things in bounds
          if (MFg2_Vector[g]-2 == MFg2_Vector[g-1]) {  // Sanity check
            gsl_complex tmp = rho_gg[g-1][gp];
            // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
            gsl_complex_mul_real(tmp, cMins_G[g]);
            left = gsl_complex_add(left, tmp);
          }
        }
        left = gsl_complex_mul_real(left, gFactor_G);
        // Is there a state with Mf` = Mf+1 and F` = F?
        if (gp + 1 < numGStates) {      // Keep things in bounds
          if (MFg2_Vector[gp]+2 == MFg2_Vector[gp+1]) {  // Sanity check
            gsl_complex tmp = rho_gg[g][gp+1];
            // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
            gsl_complex_mul_real(tmp, cPlus_G[gp]);
            rigt = gsl_complex_sub(rigt, tmp);
          }
        }
        // Is there a state with Mf` = Mf-1 and F` = F?
        if (gp > 0) {                   // Keep things in bounds
          if (MFg2_Vector[gp]-2 == MFg2_Vector[gp-1]) {  // Sanity check
            gsl_complex tmp = rho_gg[g][gp-1];
            // printf("%g + %g i\n", GSL_REAL(tmp), GSL_IMAG(tmp));
            gsl_complex_mul_real(tmp, cMins_G[gp]);
            rigt = gsl_complex_sub(rigt, tmp);
          }
        }
        rigt = gsl_complex_mul_real(rigt, gFactor_G);
        Bx = gsl_complex_add(left, rigt);
        Bx = gsl_complex_mul_imag(Bx,
                                  -_bohr_magneton*eigen.field.B_x
                                  /2.0/_planck_hbar);
        // printf("Bx contribution is %g + %g i\n", GSL_REAL(Bx), GSL_IMAG(Bx));
        drho_gg[g][gp] = gsl_complex_add(drho_gg[g][gp], Bx);
      }  // End if coherences
    }    // End gp
  }      // End g
  // F-Ground states rho_ff:  Equation 33
  // Zeeman coherences are when f != f`
  for (int f = 0; f < numFStates; f++) {
    for (int fp = 0; fp < numFStates; fp++) {
      gsl_complex oCoherences = gsl_complex_rect(0.0, 0.0);
      if ((f == fp) || (gs_Zeeman)) {
        for (int q = 0; q < 3; q += 2) {  // Only q = 0, 2 will have lasers
          gsl_complex qTerm = gsl_complex_rect(0.0, 0.0);
          for (int epp = 0; epp < numEStates; epp++) {
            gsl_complex left, right;  // First and second halves in
            // top line of equation 34
            left = gsl_complex_rect(pow(-1.0, (0))*
                                    dipole_moment_ef[epp][f]*a_ef[epp][f][q],
                                    0.0);
            // The 2-q is added due to de Clerq equation I.2
            left = gsl_complex_mul(left,
                                   delta_ef[epp][fp]);
            right = gsl_complex_rect(dipole_moment_ef[epp][fp]*a_ef[epp][fp][q],
                                     0.0);
            right = gsl_complex_mul(right,
                                    gsl_complex_conjugate(delta_ef[epp][f]));
            left = gsl_complex_sub(left, right);
            qTerm = gsl_complex_add(qTerm, left);
          }  // End epp-loop
          qTerm = gsl_complex_mul_real(qTerm, laser_fe.field[q]);
          oCoherences = gsl_complex_add(oCoherences, qTerm);
        }  // End q-loop
        oCoherences = gsl_complex_mul_imag(oCoherences,
                                           1.0/(2*_planck_hbar));
        drho_ff[f][fp] = gsl_complex_add(drho_ff[f][fp], oCoherences);
        { // Spontaneous decay term
          gsl_complex spon_decay = gsl_complex_rect(0.0, 0.0);
          for (int e = 0; e < numEStates; e++) {
            for (int ep = 0; ep < numEStates; ep++) {
              if (MFe2_Vector[e] - MFe2_Vector[ep] ==
                  MFf2_Vector[f] - MFf2_Vector[fp]) {
                int q = MFe2_Vector[e] - MFf2_Vector[f];
                if (q%2 != 0) {
                  printf("F - spontaneous decay: q2 not even\n");
                  exit(1);
                }
                q = (q/2) + 1;    // Converts q from twice the difference in M_f
                // to an appropriate index
                double coupling = a_ef[e][f][q] * a_ef[ep][fp][q];
                gsl_complex left = gsl_complex_mul_real(rho_ee[e][ep],
                                                        coupling);
                spon_decay = gsl_complex_add(spon_decay, left);
              }  // End if M_f sublevels satisfy m_e - m_e` = m_f - m_f`
            }    // End ep sum
          }      // End e sum
          spon_decay = gsl_complex_mul_real(spon_decay, ((1.0/data.atom.tau)));
          double angFrequency = 2*M_PI*(nu_F[fp] - nu_F[f]);
          gsl_complex right = gsl_complex_mul_imag(rho_ff[f][fp],
                                                   angFrequency);
          spon_decay = gsl_complex_sub(spon_decay, right);
          drho_ff[f][fp] = gsl_complex_add(drho_ff[f][fp], spon_decay);
        }  // End spontaneous decay term
        // Get ready for transverse magnetic field!
        gsl_complex Bx = gsl_complex_rect(0.0, 0.0);
        gsl_complex left = gsl_complex_rect(0.0, 0.0);
        gsl_complex rigt = gsl_complex_rect(0.0, 0.0);
        // printf("For e = %d, comparing %d to %d\n", e, MFe2_Vector[e],
        //        MFe2_Vector[e+1]);
        // Is there a state with Mf` = Mf+1 and F` = F?
        if (f+1 < numFStates) {         // Keep things in bounds
          if (MFf2_Vector[f]+2 == MFf2_Vector[f+1]) {  // Sanity check
            gsl_complex tmp = rho_ff[f+1][fp];
            // printf("Bx = %g + %g i + ", GSL_REAL(tmp), GSL_IMAG(tmp));
            gsl_complex_mul_real(tmp, cPlus_F[f]);
            left = gsl_complex_add(left, tmp);
          }
        }
        // Is there a state with Mf` = Mf-1 and F` = F?
        if (f > 0) {                    // Keep things in bounds
          if (MFf2_Vector[f]-2 == MFf2_Vector[f-1]) {  // Sanity check
            gsl_complex tmp = rho_ff[f-1][fp];
            // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
            gsl_complex_mul_real(tmp, cMins_F[f]);
            left = gsl_complex_add(left, tmp);
          }
        }
        left = gsl_complex_mul_real(left, gFactor_F);
        // Is there a state with Mf` = Mf+1 and F` = F?
        if (fp+1 < numFStates) {        // Keep things in bounds
          if (MFf2_Vector[fp]+2 == MFf2_Vector[fp+1]) {  // Sanity check
            gsl_complex tmp = rho_ff[f][fp+1];
            // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
            gsl_complex_mul_real(tmp, cPlus_F[fp]);
            rigt = gsl_complex_sub(rigt, tmp);
          }
        }
        // Is there a state with Mf` = Mf-1 and F` = F?
        if (fp > 0) {                   // Keep things in bounds
          if (MFf2_Vector[fp]-2 == MFf2_Vector[fp-1]) {  // Sanity check
            gsl_complex tmp = rho_ff[f][fp-1];
            // printf("%g + %g i\n", GSL_REAL(tmp), GSL_IMAG(tmp));
            gsl_complex_mul_real(tmp, cMins_F[fp]);
            rigt = gsl_complex_sub(rigt, tmp);
          }
        }
        rigt = gsl_complex_mul_real(rigt, gFactor_F);
        Bx = gsl_complex_add(left, rigt);
        Bx = gsl_complex_mul_imag(Bx,
                                  -_bohr_magneton*eigen.field.B_x
                                  /2.0/_planck_hbar);
        // printf("Bx contribution is %g + %g i\n", GSL_REAL(Bx), GSL_IMAG(Bx));
        drho_ff[f][fp] = gsl_complex_add(drho_ff[f][fp], Bx);
      }  // End if doing coherences
    }    // End fp loop
  }      // End f loop
  // EG-Optical Coherences delta_eg:  Equation 36
  bool debug_eg = false;
  for (int e = 0; e < numEStates; e++) {
    for (int g = 0; g < numGStates; g++) {
      if (debug_eg) printf("e = %d, g = %d\n", e, g);
      gsl_complex gLaser_term = gsl_complex_rect(0.0, 0.0);
      for (int q = 0; q < 3; q += 2) {  // Only q = 0, 2 will have lasers
          if (debug_eg) printf("\t q = %d\n", q);
          gsl_complex gsum = gsl_complex_rect(0.0, 0.0);
          gsl_complex esum = gsl_complex_rect(0.0, 0.0);
          for (int gpp = 0; gpp < numGStates; gpp++) {
            gsl_complex temp = gsl_complex_rect(dipole_moment_eg[e][gpp], 0.0);
            temp = gsl_complex_mul_real(temp, a_eg[e][gpp][q]);
            temp = gsl_complex_mul(temp, rho_gg[gpp][g]);
            gsum = gsl_complex_add(gsum, temp);
          }  // End gpp-loop
          for (int epp = 0; epp < numEStates; epp++) {
            gsl_complex temp = gsl_complex_rect(dipole_moment_eg[epp][g], 0.0);
            temp = gsl_complex_mul_real(temp, a_eg[epp][g][q]);
            temp = gsl_complex_mul(temp, rho_ee[e][epp]);
            esum = gsl_complex_add(esum, temp);
          }  // End epp-loop
          gsum = gsl_complex_sub(gsum, esum);
          gsum = gsl_complex_mul_real(gsum, laser_ge.field[q]);
          gLaser_term = gsl_complex_add(gLaser_term, gsum);
      }  // End q-loop
      gLaser_term = gsl_complex_mul_imag(gLaser_term, 1.0/(2.0*_planck_hbar));
      ddelta_eg[e][g] = gsl_complex_add(ddelta_eg[e][g], gLaser_term);

      // GS-Hyperfine coherence term
      if (gs_hyperfine) {
        gsl_complex dotProduct = gsl_complex_rect(0.0, 0.0);
        for (int q = 0; q < 3; q+=2) {  // Only the q = 0, 2 will have lasers
          gsl_complex fsum = gsl_complex_rect(0.0, 0.0);
          for (int fpp = 0; fpp < numFStates; fpp++) {
            gsl_complex temp = gsl_complex_rect(dipole_moment_ef[e][fpp], 0.0);
            temp = gsl_complex_mul(temp, gsl_complex_rect(a_ef[e][fpp][q],
                                                          0.0));
            temp = gsl_complex_mul(temp, delta_fg[fpp][g]);
            fsum = gsl_complex_add(fsum, temp);
          }  // End loop over gpp
          fsum = gsl_complex_mul(fsum,
                                 gsl_complex_rect(laser_fe.field[q], 0.0));
          dotProduct = gsl_complex_add(dotProduct, fsum);
        }   // End loop over q
        gsl_complex hyperfine_term = gsl_complex_mul_imag(dotProduct,
                                                          0.5/_planck_hbar);
        ddelta_eg[e][g] = gsl_complex_add(ddelta_eg[e][g], hyperfine_term);
      }  // End if gs_hyperfine

      // Detuning term -i(\omega_eg - \omega_1)*\delta_eg
      double detune_d = -2*M_PI*((nu_E[e] - nu_G[g]) - laser_ge.nu);
      double linewidth = -0.5*((1.0/data.atom.tau) +
                               (2*M_PI*laser_ge.linewidth));
      gsl_complex detune_g = gsl_complex_rect(linewidth, detune_d);
      detune_g = gsl_complex_mul(detune_g, delta_eg[e][g]);
      ddelta_eg[e][g] = gsl_complex_add(ddelta_eg[e][g], detune_g);

      // Get ready for transverse magnetic field!
      gsl_complex Bx = gsl_complex_rect(0.0, 0.0);
      gsl_complex left = gsl_complex_rect(0.0, 0.0);
      gsl_complex rigt = gsl_complex_rect(0.0, 0.0);
      // printf("For e = %d, comparing %d to %d\n", e, MFe2_Vector[e],
      //        MFe2_Vector[e+1]);
      // Is there a state with Mf` = Mf+1 and F` = F?
      if (e+1 < numEStates) {           // Keep things in bounds
        if (MFe2_Vector[e]+2 == MFe2_Vector[e+1]) {  // Sanity check
          gsl_complex tmp = delta_eg[e+1][g];
          // printf("Bx = %g + %g i + ", GSL_REAL(tmp), GSL_IMAG(tmp));
          gsl_complex_mul_real(tmp, cPlus_E[e]);
          left = gsl_complex_add(left, tmp);
        }
      }
      // Is there a state with Mf` = Mf-1 and F` = F?
      if (e > 0) {                      // Keep things in bounds
        if (MFe2_Vector[e]-2 == MFe2_Vector[e-1]) {  // Sanity check
          gsl_complex tmp = delta_eg[e-1][g];
          // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
          gsl_complex_mul_real(tmp, cMins_E[e]);
          left = gsl_complex_add(left, tmp);
        }
      }
      left = gsl_complex_mul_real(left, gFactor_E[e]);
      // Is there a state with Mf` = Mf+1 and F` = F?
      if (g+1 < numGStates) {           // Keep things in bounds
        if (MFg2_Vector[g]+2 == MFg2_Vector[g+1]) {  // Sanity check
          gsl_complex tmp = delta_eg[e][g+1];
          // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
          gsl_complex_mul_real(tmp, cPlus_G[g]);
          rigt = gsl_complex_sub(rigt, tmp);
        }
      }
      // Is there a state with Mf` = Mf-1 and F` = F?
      if (g > 0) {                      // Keep things in bounds
        if (MFg2_Vector[g]-2 == MFg2_Vector[g-1]) {  // Sanity check
          gsl_complex tmp = delta_eg[e][g-1];
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
      ddelta_eg[e][g] = gsl_complex_add(ddelta_eg[e][g], Bx);
    }  //  End g
  }    // End e

  // EF-Optical Coherences delta_ef:  Equation 35
  bool debug_ef = false;
  for (int e = 0; e < numEStates; e++) {
    for (int f = 0; f < numFStates; f++) {
      if (debug_ef) printf("e = %d, f = %d\n", e, f);
      gsl_complex fLaser_term = gsl_complex_rect(0.0, 0.0);
      for (int q = 0; q < 3; q +=2) {  // Only q = 0, 2 will have lasers
          if (debug_ef) printf("\t q = %d\n", q);
          gsl_complex fsum = gsl_complex_rect(0.0, 0.0);
          gsl_complex esum = gsl_complex_rect(0.0, 0.0);
          for (int fpp = 0; fpp < numFStates; fpp++) {
            gsl_complex temp = gsl_complex_rect(dipole_moment_ef[e][fpp], 0.0);
            temp = gsl_complex_mul_real(temp, a_ef[e][fpp][q]);
            temp = gsl_complex_mul(temp, rho_ff[fpp][f]);
            fsum = gsl_complex_add(fsum, temp);
          }  // End fpp-loop
          for (int epp = 0; epp < numEStates; epp++) {
            gsl_complex temp = gsl_complex_rect(dipole_moment_ef[epp][f], 0.0);
            temp = gsl_complex_mul_real(temp, a_ef[epp][f][q]);
            temp = gsl_complex_mul(temp, rho_ee[e][epp]);
            esum = gsl_complex_add(esum, temp);
          }  // End epp-loop
          fsum = gsl_complex_sub(fsum, esum);
          fsum = gsl_complex_mul_real(fsum, laser_fe.field[q]);
          fLaser_term = gsl_complex_add(fLaser_term, fsum);
      }  // End q-loop
      fLaser_term = gsl_complex_mul_imag(fLaser_term, 1.0/(2.0*_planck_hbar));
      ddelta_ef[e][f] = gsl_complex_add(ddelta_ef[e][f], fLaser_term);

      // GS-Hyperfine coherence term
      if (gs_hyperfine) {
        gsl_complex dotProduct = gsl_complex_rect(0.0, 0.0);
        for (int q = 0; q < 3; q+=2) {  // Only the q = 0, 2 will have lasers
          gsl_complex gsum = gsl_complex_rect(0.0, 0.0);
          for (int gpp = 0; gpp < numGStates; gpp++) {
            gsl_complex temp = gsl_complex_rect(dipole_moment_eg[e][gpp], 0.0);
            temp = gsl_complex_mul(temp, gsl_complex_rect(a_eg[e][gpp][q],
                                                          0.0));
            temp = gsl_complex_mul(temp,
                                  gsl_complex_conjugate(delta_fg[f][gpp]));
            gsum = gsl_complex_add(gsum, temp);
          }  // End loop over gpp
          gsum = gsl_complex_mul(gsum,
                                 gsl_complex_rect(laser_ge.field[q], 0.0));
          dotProduct = gsl_complex_add(dotProduct, gsum);
        }   // End loop over q
        gsl_complex hyperfine_term = gsl_complex_mul_imag(dotProduct,
                                                          0.5/_planck_hbar);
        ddelta_ef[e][f] = gsl_complex_add(ddelta_ef[e][f], hyperfine_term);
      }  // End if gs_hyperfine

      // Detuning term -i(\omega_ef - \omega_1)*\delta_ef
      double detune_d = -2*M_PI*((nu_E[e] - nu_F[f]) - laser_fe.nu);
      double linewidth = -0.5*((1.0/data.atom.tau) +
                               (2.0*M_PI*laser_fe.linewidth));
      gsl_complex detune_f = gsl_complex_rect(linewidth, detune_d);
      detune_f = gsl_complex_mul(detune_f, delta_ef[e][f]);
      ddelta_ef[e][f] = gsl_complex_add(ddelta_ef[e][f], detune_f);

      // Get ready for transverse magnetic field!
      gsl_complex Bx = gsl_complex_rect(0.0, 0.0);
      gsl_complex left = gsl_complex_rect(0.0, 0.0);
      gsl_complex rigt = gsl_complex_rect(0.0, 0.0);
      // printf("For e = %d, comparing %d to %d\n", e, MFe2_Vector[e],
      //        MFe2_Vector[e+1]);
      // Is there a state with Mf` = Mf+1 and F` = F?
      if (e+1 < numEStates) {           // Keep things in bounds
        if (MFe2_Vector[e]+2 == MFe2_Vector[e+1]) {  // Sanity check
          gsl_complex tmp = delta_ef[e+1][f];
          // printf("Bx = %g + %g i + ", GSL_REAL(tmp), GSL_IMAG(tmp));
          gsl_complex_mul_real(tmp, cPlus_E[e]);
          left = gsl_complex_add(left, tmp);
        }
      }
      // Is there a state with Mf` = Mf-1 and F` = F?
      if (e > 0) {                      // Keep things in bounds
        if (MFe2_Vector[e]-2 == MFe2_Vector[e-1]) {  // Sanity check
          gsl_complex tmp = delta_ef[e-1][f];
          // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
          gsl_complex_mul_real(tmp, cMins_E[e]);
          left = gsl_complex_add(left, tmp);
        }
      }
      left = gsl_complex_mul_real(left, gFactor_E[e]);
      // Is there a state with Mf` = Mf+1 and F` = F?
      if (f+1 < numFStates) {           // Keep things in bounds
        if (MFf2_Vector[f]+2 == MFf2_Vector[f+1]) {  // Sanity check
          gsl_complex tmp = delta_ef[e][f+1];
          // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
          gsl_complex_mul_real(tmp, cPlus_F[f]);
          rigt = gsl_complex_sub(rigt, tmp);
        }
      }
      // Is there a state with Mf` = Mf-1 and F` = F?
      if (f > 0) {                      // Keep things in bounds
        if (MFf2_Vector[f]-2 == MFf2_Vector[f-1]) {  // Sanity check
          gsl_complex tmp = delta_ef[e][f-1];
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
      ddelta_ef[e][f] = gsl_complex_add(ddelta_ef[e][f], Bx);
    }  //  End f
  }    // End e

  // FG-Ground state hyperfine coherences: Equation 37
  if (gs_hyperfine) {
    bool debug_fg = false;
    for (int f = 0; f < numFStates; f++) {
      for (int g = 0; g < numGStates; g++) {
        if (debug_fg) printf("f = %d, g = %d\n", f, g);
        ddelta_fg[f][g] = gsl_complex_rect(0.0, 0.0);
        // F-laser term & G-Laser term
        gsl_complex fLaser_term = gsl_complex_rect(0.0, 0.0);
        {                                // First sum in equation 37
          for (int q = 0; q < 3; q += 2) {  // Only q = 0, 2 3ill have lasers
            if (debug_fg) printf("\t q = %d\n", q);
            gsl_complex esum = gsl_complex_rect(0.0, 0.0);
            for (int epp = 0; epp < numEStates; epp++) {
              if (debug_fg) printf("\t\t epp = %d\n", epp);
              gsl_complex temp = gsl_complex_rect(dipole_moment_ef[epp][f]*
                                                  pow(-1.0, 0.0), 0.0);
              temp = gsl_complex_mul(temp, gsl_complex_rect(a_ef[epp][f][q],
                                                            0.0));
              temp = gsl_complex_mul(temp, delta_eg[epp][g]);
              esum = gsl_complex_add(esum, temp);
              if (debug_fg) {
                printf("\t\t\t  temp = %8.6G + %8.6G i", GSL_REAL(temp),
                       GSL_IMAG(temp));
                printf(" ---> esum = %8.6G + %8.6G i\n", GSL_REAL(esum),
                       GSL_IMAG(esum));
              }  // End debug
            }  // End loop of epp (f-laser term)
            esum = gsl_complex_mul(esum,
                                   gsl_complex_rect(laser_fe.field[q], 0.0));
            fLaser_term = gsl_complex_add(fLaser_term, esum);
            if (debug_fg) {
              printf("\t\t fLaser_term = %8.6G + %8.6G i\n",
                     GSL_REAL(fLaser_term), GSL_IMAG(fLaser_term));
            }
          }  // End loop over q = 0, 2
        }   // End first sum in equation 37
        gsl_complex gLaser_term = gsl_complex_rect(0.0, 0.0);
        {                                // Second sum in equation 37
          for (int q = 0; q < 3; q += 2) {  // Only q = 0, 2 will have lasers
            if (debug_fg) printf("\t q = %d\n", q);
            gsl_complex esum = gsl_complex_rect(0.0, 0.0);
            for (int epp = 0; epp < numEStates; epp++) {
              if (debug_fg) printf("\t\t epp = %d\n", epp);
              gsl_complex temp = gsl_complex_rect(dipole_moment_eg[epp][g],
                                                  0.0);
              temp = gsl_complex_mul(temp, gsl_complex_rect(a_eg[epp][g][q],
                                                            0.0));
              temp = gsl_complex_mul(temp,
                                     gsl_complex_conjugate(delta_ef[epp][f]));
              esum = gsl_complex_add(esum, temp);
              if (debug_fg) {
                printf("\t\t\t  temp = %8.6G + %8.6G i", GSL_REAL(temp),
                       GSL_IMAG(temp));
                printf(" ---> esum = %8.6G + %8.6G i\n", GSL_REAL(esum),
                       GSL_IMAG(esum));
              }  // End debug
            }  // End loop of epp (f-laser term)
            esum = gsl_complex_mul(esum,
                                   gsl_complex_rect(laser_ge.field[q], 0.0));
            gLaser_term = gsl_complex_add(gLaser_term, esum);
            if (debug_fg) {
              printf("\t\t gLaser_term = %8.6G + %8.6G i\n",
                     GSL_REAL(gLaser_term), GSL_IMAG(gLaser_term));
            }
          }  // End loop over q = 0, 2
        }    // End first sum in equation 37

        // Now the `f-Laser' term will stand for the subtraction of f-g
        fLaser_term = gsl_complex_sub(fLaser_term, gLaser_term);
        fLaser_term = gsl_complex_mul_imag(fLaser_term, 0.5/_planck_hbar);
        if (debug_fg) {
          printf("Total laser term = %8.6G + %8.6G i\n", GSL_REAL(fLaser_term),
                 GSL_IMAG(fLaser_term));
        }
        ddelta_fg[f][g] = gsl_complex_add(ddelta_fg[f][g], fLaser_term);
        // Now the line width and angular frequencies terms
        // Here, T&J write that the ground state hyperfine coherences decay
        // with a rate 0.5*(gamma_1+gamma_2).  Note that this is for
        // uncorrelated laser beams.  In typical TRINAT application, the beams
        // are 100% correlated and this term drops out.  More accurately, T&J
        // leave off the correlation term -(2*gamma_12) between the two lasers
        // as per their assumptions.  Therefore I am leaving this term out of
        // the calculation.  Note that all other terms remain
        // Useful references for this point:
        // S. Gu et. al. Op. Comm.  220 (2003) 365-370
        // E. Arimondo, Progress in Optics XXXV (1996) 257 (eq 2.37)
        // B. Dalton and P. Knight J. Phys. B 15 (1982) 3997-4015 (fig 2)
        // double real = M_PI*(laser_fe.linewidth + laser_ge.linewidth);
        double real = 500.0 * _Hz;
        double imag = 2*M_PI*((fabs(nu_F[f] - nu_G[g]))
                              - (laser_ge.nu - laser_fe.nu));
        gsl_complex decayTerm = gsl_complex_rect(real, imag);
        decayTerm = gsl_complex_mul(decayTerm, delta_fg[f][g]);
        if (debug_fg) {
          printf("Decay term = %8.6G + %8.6G i\n", GSL_REAL(decayTerm),
                 GSL_IMAG(decayTerm));
        }
        ddelta_fg[f][g] = gsl_complex_sub(ddelta_fg[f][g], decayTerm);

      // Get ready for transverse magnetic field!
      gsl_complex Bx = gsl_complex_rect(0.0, 0.0);
      gsl_complex left = gsl_complex_rect(0.0, 0.0);
      gsl_complex rigt = gsl_complex_rect(0.0, 0.0);
      // printf("For e = %d, comparing %d to %d\n", e, MFe2_Vector[e],
      //        MFe2_Vector[e+1]);
      // Is there a state with Mf` = Mf+1 and F` = F?
      // if (f+1 < numFStates) {           // Keep things in bounds
      //   if (MFf2_Vector[f]+2 == MFf2_Vector[f+1]) {  // Sanity check
      //     gsl_complex tmp = delta_fg[f+1][g];
      //     // printf("Bx = %g + %g i + ", GSL_REAL(tmp), GSL_IMAG(tmp));
      //     gsl_complex_mul_real(tmp, cPlus_F[f]);
      //     left = gsl_complex_add(left, tmp);
      //   }
      // }
      // // Is there a state with Mf` = Mf-1 and F` = F?
      // if (f > 0) {                      // Keep things in bounds
      //   if (MFf2_Vector[f]-2 == MFf2_Vector[f-1]) {  // Sanity check
      //     gsl_complex tmp = delta_fg[f-1][g];
      //     // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
      //     gsl_complex_mul_real(tmp, cMins_F[f]);
      //     left = gsl_complex_add(left, tmp);
      //   }
      // }
      // left = gsl_complex_mul_real(left, gFactor_F);
      // // Is there a state with Mf` = Mf+1 and F` = F?
      // if (g+1 < numGStates) {           // Keep things in bounds
      //   if (MFg2_Vector[g]+2 == MFg2_Vector[g+1]) {  // Sanity check
      //     gsl_complex tmp = delta_fg[f][g+1];
      //     // printf("%g + %g i - ", GSL_REAL(tmp), GSL_IMAG(tmp));
      //     gsl_complex_mul_real(tmp, cPlus_G[g]);
      //     rigt = gsl_complex_sub(rigt, tmp);
      //   }
      // }
      // // Is there a state with Mf` = Mf-1 and F` = F?
      // if (g > 0) {                      // Keep things in bounds
      //   if (MFg2_Vector[g]-2 == MFg2_Vector[g-1]) {  // Sanity check
      //     gsl_complex tmp = delta_fg[f][g-1];
      //     // printf("%g + %g i\n", GSL_REAL(tmp), GSL_IMAG(tmp));
      //     gsl_complex_mul_real(tmp, cMins_G[g]);
      //     rigt = gsl_complex_sub(rigt, tmp);
      //   }
      // }
      // rigt = gsl_complex_mul_real(rigt, gFactor_G);
      Bx = gsl_complex_add(left, rigt);
      Bx = gsl_complex_mul_imag(Bx,
                                -_bohr_magneton*eigen.field.B_x
                                /2.0/_planck_hbar);
      // printf("Bx contribution is %g + %g i\n", GSL_REAL(Bx), GSL_IMAG(Bx));
      // ddelta_fg[f][g] = gsl_complex_add(ddelta_fg[f][g], Bx);
      }   // End loop over g states
    }     // End loop over f states
  }       // End if doing gs_hyperine coherences

  for (int e = 0; e < numEStates; e++) {
    for (int ep = 0; ep < numEStates; ep++) {
      drho_ee[e][ep] = gsl_complex_mul_real(drho_ee[e][ep], dt);
      rho_ee[e][ep] = gsl_complex_add(rho_ee[e][ep], drho_ee[e][ep]);
      if ((e == ep) && GSL_REAL(rho_ee[e][ep]) < 0.0) {
        rho_ee[e][ep] = gsl_complex_rect(0.0, 0.0);
      }
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
      if ((f == fp) && GSL_REAL(rho_ff[f][fp]) < 0.0) {
        rho_ff[f][fp] = gsl_complex_rect(0.0, 0.0);
      }
    }
    for (int g = 0; g < numGStates; g++) {
      ddelta_fg[f][g] = gsl_complex_mul_real(ddelta_fg[f][g], dt);
      delta_fg[f][g] = gsl_complex_add(delta_fg[f][g], ddelta_fg[f][g]);
    }
  }
  for (int g = 0; g < numGStates; g++) {
    for (int gp = 0; gp < numGStates; gp++) {
      drho_gg[g][gp] = gsl_complex_mul_real(drho_gg[g][gp], dt);
      rho_gg[g][gp] = gsl_complex_add(rho_gg[g][gp], drho_gg[g][gp]);
      if ((g == gp) && GSL_REAL(rho_gg[g][gp]) < 0.0) {
        rho_gg[g][gp] = gsl_complex_rect(0.0, 0.0);
      }
    }
  }
}

void Density_Matrix::setup_dipole_moments(double gamma) {
  for (int e = 0; e < numEStates; e++) {
    double ang_freq_ex = nu_E[e] * 2.0 * M_PI;
    for (int g = 0; g < numGStates; g++) {
      if (op_verbose) printf("e = %d <---> g = %d \t", e, g);
      double ang_freq_gr = nu_G[g] * 2.0 * M_PI;
      dipole_moment_eg[e][g] = set_dipole_moment(gamma, ang_freq_ex,
                                                ang_freq_gr);
    }
    for (int f = 0; f < numFStates; f++) {
      if (op_verbose) printf("e = %d <---> f = %d\t", e, f);
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
  if (op_verbose) {
    printf("gamma = %5.3G MHz\t delta_omega = %19.16G MHz\t", gamma/_MHz,
           fabs(omega_ex-omega_gr)/_MHz);
    printf("dipole moment = %8.6G e*nm\t", dipole/(_elementary_charge*_nm));
    printf("laser rate = %8.6G MHz\n",
           (dipole*data.laser_fe.field[2]/(2*_planck_hbar))/_MHz);
  }
  return dipole;
}

void Density_Matrix::print_rabi_frequencies(FILE * des) {
  // Calculate Rabi frequencies based on Metcalf equation 1.10.  Note that this
  // equation has the angular frequency!  I will also print the period to avoid
  // confusion.
  fprintf(des, "*****Rabi Frequencies*****\n");
  for (int q = 0; q < 3; q++) {
    if (laser_fe.field[q] > 0.0 || laser_ge.field[q] > 0.0)
      fprintf(des, "***** q = %d *****\n", q-1);
    for (int e = 0; e < numEStates; e++) {
      if (laser_ge.field[q] > 0.0) {
        for (int g = 0; g < numGStates; g++) {
          double rabi = dipole_moment_eg[e][g] * a_eg[e][g][q] *
            laser_ge.field[q];
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
      if (laser_fe.field[q] > 0.0) {
        for (int f = 0; f < numFStates; f++) {
          double rabi = dipole_moment_ef[e][f] * a_ef[e][f][q] *
            laser_fe.field[q];
          rabi = rabi / _planck_hbar;
          double detune = 2*M_PI*(laser_fe.nu - (nu_E[e] - nu_F[f]));
          rabi = sqrt(pow(rabi, 2.0) + pow(detune, 2.0));
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

void Density_Matrix::change_magnetic_field(double newfield) {
  OpticalPumping_Method::change_magnetic_field(newfield);
  setup_dipole_moments(1.0/eigen.atom.tau);
}
