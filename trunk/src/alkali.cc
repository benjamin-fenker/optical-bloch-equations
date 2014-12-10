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

  if (strcmp(isotope.c_str(), "K37") == 0 ||
      strcmp(isotope.c_str(), "37K") == 0 ) {
    // See Besch 1968 and DM Thesis
    set_I2 = 3;
    set_Aj_g = 120.135 *_MHz;           // Williamson97
    set_Aj_e = 14.3 *_MHz;              // Williamson97
    set_g_I = .2029;                    // vonPlaten71
    // set_tau = (26.2 * _ns);          // From Melconian05 for P3/2
    set_tau = (1./0.382E+08)*_s;        // 26.2 ns from JB's code
    // double lambda = 769.9 * _nm;     // From Melconian05
    double lambda = (1./12985.) * _cm;  // 770 nm from JB's code

    if (Je2 == 3) {
      lambda = 766.7 * _nm;  // Melconian05
      set_Aj_e = 8.05 *_MHz; // Williamson97
    }

    set_nu_excited = (_speed_of_light / lambda);

  } else if (strcmp(isotope.c_str(), "K38") == 0||
             strcmp(isotope.c_str(), "38K") == 0 ) {
    // See Besch 1968 and DM Thesis
    set_I2 = 6;
    set_Aj_g = 707.65 *_MHz;            // Williamson97
    set_Aj_e = 85.3 *_MHz;              // Williamson97
    set_g_I = .2029;                    // Copied from 37K
    set_tau = (1./0.382E+08)*_s;        // 26.2 ns from JB's code
    double lambda = 770.1098 * _nm;     // Williamson97

    if (Je2 == 3) {
      lambda = 766.7017 * _nm;          // Williamson97
      set_Aj_e = 55.75 *_MHz;           // Williamson97
    }

    set_nu_excited = (_speed_of_light / lambda);
  } else if (strcmp(isotope.c_str(), "K47") == 0||
             strcmp(isotope.c_str(), "47K") == 0 ) {
    // MAKE UP VALUES, JUST NEEDED AN ISOTOPE FOR TESTING WITH I = 1/2
    set_I2 = 1;
    set_Aj_g = 707.65*_MHz;
    set_Aj_e = 85.3*_MHz;
    set_g_I = .2029;
    set_tau = (1./0.382E+08)*_s;        // From JB's code 26.2 ns
    double lambda = 769.9 * _nm;
    if (Je2 == 3) {
      lambda = 766.7 * _nm;
      set_Aj_e = 8.05 *_MHz;
    }
    set_nu_excited = (_speed_of_light / lambda);

  } else if (strcmp(isotope.c_str(), "K41") == 0 ||
            strcmp(isotope.c_str(), "41K") == 0 ) {
    set_I2 = 3;                 // Nuclear spin 3/2
    set_Aj_g = (254.01/2.0)*_MHz;       // Tiecke10
    set_Aj_e = (30.5 /2.0)*_MHz;        // Tiecke10
    set_g_I = -0.00007790600;           // Tiecke10
    set_tau = 26.7*_ns;                 // Wang1997
    //set_tau = (1./0.382E+08)*_s;  // JB's code
    double lambda = 770.107919 *_nm;    // Tiecke10
    //double lambda = (1./12985.) * _cm;  // 770 nm from JB's code
    if (Je2 == 3) {
      set_tau = 26.4 *_ns;               // Tiecke10
      lambda = 766.7004587 *_ns;        // Tiecke10
      set_Aj_e = 8.4 * _MHz;            // Tiecke10
    }
    set_nu_excited = (_speed_of_light / lambda);
  } else if (strcmp(isotope.c_str(), "Rb85") == 0 ||
             strcmp(isotope.c_str(), "85Rb") == 0) {
    set_I2 = 5;
    set_Aj_g = 1011.9108130*_MHz;       // Arimondo1977
    set_Aj_e = 120.640*_MHz;            // Banerjee2004
    set_g_I = -0.00029364000;           // Arimondo1977
    set_tau = 27.75*_ns;                // Gutterres2002
    double lambda = 794.979014933*_nm;  // Steck2008
    if (Je2 == 3) {
      set_Aj_e = 25.009*_MHz;           // Arimondo1977
      lambda = 780.24136827127*_nm;     // Steck2008
      set_tau = 26.25*_ns;              // Gutterres2002
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

double Alkali::getLaserFrequency(atom_data atom, magnetic_field_data field,
                                 int Fg2, int Mfg2, int Fe2, int Mfe2,
                                 double detune) {
  bool debug = false;
  OpticalPumping_Method *opm = new Rate_Equations();
  //  Eigenvector_Helper *eigen = new Eigenvector_Helper(atom, field);
  double energy_gr = opm -> set_frequency(0.0, atom.I2, 1, Fg2, Mfg2, 0,
                                       atom.Aj_g, atom.g_I, field.B_z);
  if (debug) printf("Ground energy: %8.6f\n", energy_gr);
  if (debug) printf("Aj_e = %g MHz\n", atom.Aj_e/_MHz);
  double energy_ex = opm -> set_frequency(atom.nu_excited, atom.I2, atom.Je2,
                                       Fe2, Mfe2, 2, atom.Aj_e, atom.g_I,
                                       field.B_z);
  if (debug) printf("Excited energy: %8.6f\n", energy_ex);
  if (debug) printf("Energy difference [MHz] = %8.6f\n",
                    (energy_ex - energy_gr)/_MHz);
  delete opm;
  return energy_ex - energy_gr + detune;
}

double Alkali::getGamma(double tau, double laser_power) {
  // No power broadening yet.  Which I_laser do I use if they are different?
  // Which lambda do I use; they will be different!
  double gamma = (1.0 / tau);
  // double one = laser_power / laser_power;
  // gamma *= one;
  return gamma;
}
