// Copyright 2012 Benjamin Fenker

#ifndef INCLUDE_OPTICAL_PUMPING_DATA_STRUCTURES_H_
#define INCLUDE_OPTICAL_PUMPING_DATA_STRUCTURES_H_

#include <vector>

using std::vector;
// Struct to hold atomic data and make function calls clearer
struct atom_data {
  int numEStates, numFStates, numGStates;
  int I2, Je2;
  double nu_excited, tau, gamma_spon;
  double Aj_g, Aj_e, g_I;
  double linewidth;
};

struct magnetic_field_data {
  double mu_B;
  double B_z;
  double B_x;
};

struct coherence_flags {
  bool zCoherences;
  bool hfCoherences_ex;
  bool hfCoherences_gr;
};

class Laser_data {
 private:
  void set_field_components();
  void set_intensity_components();
 public:
  Laser_data();
  // Old constructor with polarizations
  // Laser_data(double nu, double power, double detune, double linewidth,
  //           double polarization[], double tau);
  Laser_data(double nu, double power, double detune, double linewidth,
	     double s3_over_s0, double tau);
  void set_saturation_intensity(double tau);
  // and set I_sat in mW/cm^2
  double nu;  // Frequency of laser.  Energy = h * nu
  double power, intensity[2], field[2];  // sigma^- , sigma^+ for each
  // power is the total laser power while intensity is the power in each
  // component
  double stokes[4];  // Stokes parameters as in Jackson
  // I only anticipate using [0] and [3] as those are the parameters that define
  // the polarization.
  double detune;
  double linewidth;  // linewidth = FWHM
  // linewidth will represent linewidth in frequency with linewidth in energy
  // = h * linewidth (as opposed to linewidth in angular frequency with
  // energy linewidth = hbar * linewidth)
  //  double polarization[3]; Don't want a seperate parameter anymore
  double saturation_intensity;
};

struct op_data_for_gsl {
  // This struct is not designed safely!!!
  // These are used in both RE and OBE codes
  int numGStates, numFStates, numEStates;
  int gStart, gEnd, fStart, fEnd, eStart, eEnd;
  atom_data atom;
  Laser_data laser_ge, laser_fe;
  vector<vector<vector<double> > > a_eg, a_ef;
  vector<double> nu_E, nu_F, nu_G;

  // These are required for RE but not allowed in OBE!!
  vector<vector<double> > transition_rate_eg, transition_rate_ef;

  // These are require for OBE but not allowed in RE!!
  vector<vector<double> > dipole_moment_eg, dipole_moment_ef;
};

#endif  // INCLUDE_OPTICAL_PUMPING_DATA_STRUCTURES_H_
