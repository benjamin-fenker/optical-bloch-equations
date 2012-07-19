// Copyright 2012 Benjamin Fenker

#ifndef NEWOBE_INCLUDE_OPTICAL_PUMPING_DATA_STRUCTURES_H_
#define NEWOBE_INCLUDE_OPTICAL_PUMPING_DATA_STRUCTURES_H_

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
 public:
  Laser_data();
  Laser_data(double nu, double power, double detune, double linewidth,
             double polarization[], double tau);
  void set_saturation_intensity(double tau);
  void set_E_field_circular();
  void set_E_field_linear();
  // and set I_sat in mW/cm^2
  double nu;  // Frequency of laser.  Energy = h * nu
  double power, field;
  double detune;
  double linewidth;  // linewidth = FWHM
  // linewidth will represent linewidth in frequency with linewidth in energy
  // = h * linewidth (as opposed to linewidth in angular frequency with
  // energy linewidth = hbar * linewidth)
  double polarization[3];
  double saturation_intensity;
};
#endif  // NEWOBE_INCLUDE_OPTICAL_PUMPING_DATA_STRUCTURES_H_
