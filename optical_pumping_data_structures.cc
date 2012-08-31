// Copyright 2012 Benjamin Fenker

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_const_mksa.h>
#include "include/optical_pumping_data_structures.h"
#include "include/units.h"

extern bool op_verbose;

Laser_data::Laser_data() {}

Laser_data::Laser_data(double set_nu, double set_power, double set_detune,
                       double set_linewidth, double set_polarization[],
                       double tau) :
  nu(set_nu), power(set_power), detune(set_detune), linewidth(set_linewidth) {
  for (int i = 0; i < 3; i++) polarization[i] = set_polarization[i];
  set_saturation_intensity(tau);
  // NEED TO UNDERSTAND WHY SETTING THIS TO LINEAR FIXES THINGS.  WHAT DOES
  // LASER INTENSITY MEAN???
  //  set_E_field_circular();
  set_E_field_linear();
}

void Laser_data::set_saturation_intensity(double tau) {
  double I_s = M_PI / (3.0*tau);
  I_s /= pow(_speed_of_light, 2.0);
  I_s *= _planck_h * pow(nu, 3);
  saturation_intensity = I_s;
  if (op_verbose) printf("I_s = %6.4G mW/cm^2\n", I_s/(_mW/_cm2));
}

void Laser_data::set_E_field_circular() {
  field = sqrt(power/(_epsilon_0 * _speed_of_light));
  if (op_verbose) printf("E field: %8.6G V/m\n", field/(_V/_m));
}

void Laser_data::set_E_field_linear() {
  field = sqrt((2.0*power)/(_epsilon_0 * _speed_of_light));
  if (op_verbose) printf("E field: %8.6G V/m\n", field/(_V/_m));
}
