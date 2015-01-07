// Authors: Benjamin Fenker 2013
// Copyright 2013 Benjamin Fenker

#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <string>

#include "optical_pumping_data_structures.h"

#ifndef INCLUDE_OPTICAL_PUMPING_H_
#define INCLUDE_OPTICAL_PUMPING_H_

using std::string;

class OpticalPumping {
 public:
  int pump(op_parameters params);
  int test();
};
#endif  // INCLUDE_OPTICAL_PUMPING_H_

