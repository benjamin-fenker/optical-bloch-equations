// Copyright 2012 Benjamin Fenker

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_const_mksa.h>
#include "include/optical_pumping_data_structures.h"
#include "include/units.h"

extern bool op_verbose;

Laser_data::Laser_data() {}

Laser_data::Laser_data(double set_nu, double set_power, double set_detune,
                       double set_linewidth, double set_s3_over_s0,
                       double tau) :
  nu(set_nu), power(set_power), detune(set_detune), linewidth(set_linewidth) {
  set_saturation_intensity(tau);
  // s1 and s2 concern the phase difference of E_+ and E_- light, but have
  // no relevance and I don't even think we measure them.  0.0 is arbitrary.
  // See Jackson (3rd edition) for definition of Stokes parameters
  stokes[0] = 2*set_power / (_epsilon_0*_speed_of_light);
  stokes[1] = 0.0;
  stokes[2] = 0.0;
  stokes[3] = stokes[0] * set_s3_over_s0;
  set_field_components();
  set_intensity_components();
}

void Laser_data::set_saturation_intensity(double tau) {
  double I_s = M_PI / (3.0*tau);
  I_s /= pow(_speed_of_light, 2.0);
  I_s *= _planck_h * pow(nu, 3);
  saturation_intensity = I_s;
  if (op_verbose) printf("I_s = %6.4G mW/cm^2\n", I_s/(_mW/_cm2));
}

void Laser_data::set_field_components() {
  // field[0] = sigma^- and field[1] = sigma^+ are the two basis vectors.
  // No capability at this time to do x-hat and y-hat basis
  if (stokes[0] <= 0.0) {
    printf("Stokes vectors not possible.  Aborting\n");
    exit(1);
  }
  field[0] = sqrt((stokes[0] - stokes[3]) / 2.0);  // l_z = -1
  field[1] = 0.0;                                  // l_z = 0
  field[2] = sqrt((stokes[0] + stokes[3]) / 2.0);  // l_z = +1
}

void Laser_data::set_intensity_components() {
  // intensity[0] = sigma^- and intensity[1] = sigma^+ are the two basis
  // vectors. No capability at this time to do x-hat and y-hat basis
  if ((field[0] < 0.0 || field[2] < 0.0) ||
      (field[0] <= 0.0 && field[2] <= 0.0)) {
    printf("Electric field not possible.  Aborting\n");
    exit(1);
  }
  for (int i = 0; i < 3; i++) {
    intensity[i] = 0.5 * _epsilon_0 * _speed_of_light * pow(field[i], 2.0);
    printf("intensity[%d] = %8.6G mW/cm2\n", i, intensity[i]/(_mW/_cm2));
  }
}
