// Authors: Benjamin Fenker 2013

// Copyright Benjamin Fenker 2013

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_const_mksa.h>

#include "include/alkali.h"
#include "include/units.h"

using std::string;

int Alkali::lookupParameters(string isotope, int Je2, int* I2, double* Aj_g,
                             double* Aj_e, double* g_I, double* nu_excited,
                             double *tau) {
  int status = 0;
  int set_I2 = 0;
  double set_Aj_g = 0;
  double set_Aj_e = 0;
  double set_g_I = 0;
  double set_nu_excited = 0;
  double set_tau = 0;

  if (strcmp(isotope.c_str(), "K37") == 0||
      strcmp(isotope.c_str(), "37K") == 0 ) {
    // See Besch 1968 and DM Thesis
    set_I2 = 3;
    set_Aj_g = 120.1336 *_MHz;
    set_Aj_e = 14.45 *_MHz;
    set_g_I = .39094;
    set_tau = 26.2 * _ns;
    double lambda = 769.9 * _nm;  // m

    if (Je2 == 3) {
      lambda = 766.7 * _nm;  // m
      set_Aj_e = 8.05 *_MHz;
    }

    set_nu_excited = (_speed_of_light / lambda);

  } else if (strcmp(isotope.c_str(), "K38") == 0||
             strcmp(isotope.c_str(), "38K") == 0 ) {
    // See Besch 1968 and DM Thesis
    set_I2 = 6;
    set_Aj_g = 707.65 *_MHz;
    set_Aj_e = 85.3 *_MHz;
    set_g_I = .2029;
    set_tau = 26.2 * _ns;
    double lambda = 769.9 * _nm;

    if (Je2 == 3) {
      lambda = 766.7 * _nm;
      set_Aj_e = 55.75 *_MHz;  // MHz
    }

    set_nu_excited = (_speed_of_light / lambda);
  } else if (strcmp(isotope.c_str(), "K47") == 0||
             strcmp(isotope.c_str(), "47K") == 0 ) {
    // MAKE UP VALUES, JUST NEEDED AN ISOTOPE FOR TESTING WITH I = 1/2
    set_I2 = 1;
    set_Aj_g = 707.65*_MHz;
    set_Aj_e = 85.3*_MHz;
    set_g_I = .2029;
    set_tau = 26.2 * _ns;
    double lambda = 769.9 * _nm;
    if (Je2 == 3) {
      lambda = 766.7 * _nm;
      set_Aj_e = 8.05 *_MHz;
    }
    set_nu_excited = (_speed_of_light / lambda);
  } else {
    status = 1;
  }

  *I2 = set_I2;
  *Aj_g = set_Aj_g;
  *Aj_e = set_Aj_e;
  *g_I = set_g_I;
  *nu_excited = set_nu_excited;
  *tau = set_tau;

  return status;
}

int Alkali::getNumberOfExcitedStates(int I2, int Je2) {
  int numEStates = (I2+1)*(Je2+1);
  return numEStates;
}

int Alkali::getNumberOfGroundStates_f(int I2) {
  int numFStates = I2+2;
  return numFStates;
}

int Alkali::getNumberOfGroundStates_g(int I2) { return I2; }

double Alkali::getLaserFrequency(double nu_excited, int Fg2, int I2,
                                 int Je2, int Fe2, double Aj_g,
                                 double Aj_e, double detune) {
  bool debug = false;
  // Zero detuning corresponds to a transition from F = Fg2, Mf = 0 to
  // F = FeMax, Mf = 0.

  if (debug) printf("Fg2 = %d, I2 = %d\n", Fg2, I2);
  // Hyperfine splitting of the ground-state energy
  // See DM Equatin 3.7 and note that 0.5Aj*(F*(F+1))... = 0.125Aj*(2F*(2F+2))
  double Jg2 = 1;
  double hfShift_ground = static_cast<double>(Fg2*(Fg2+2)) - (I2*(I2+2)) -
    (Jg2*(Jg2+2));
  hfShift_ground *= Aj_g/8.0;
  if (debug) printf("hfShift_ground = %8.6G MHz\t ", hfShift_ground/_MHz);

  // Hyperfine splitting of the excited-state energy
  double hfShift_excited = static_cast<double>(Fe2*(Fe2+2)) - (I2*(I2+2)) -
    (Je2*(Je2+2));
  hfShift_excited *= Aj_e/8.0;
  if (debug) printf("hfShift_excited = %8.6G MHz\n", hfShift_excited/_MHz);
  if (debug) printf("detune = %8.6G MHz\n", detune/_MHz);
  return nu_excited - hfShift_ground + hfShift_excited + detune;
}

double Alkali::getLaserFrequency(atom_data atom, int Fg2, int Fe2,
                                 double detune) {
  return getLaserFrequency(atom.nu_excited, Fg2, atom.I2, atom.Je2, Fe2,
                           atom.Aj_g, atom.Aj_e, detune);
}

double Alkali::getGamma(double tau, double laser_power) {
  // No power broadening yet.  Which I_laser do I use if they are different?
  // Which lambda do I use; they will be different!
  double gamma = 1.0 / tau;
  double one = laser_power / laser_power;
  gamma *= one;
  return gamma;
}
