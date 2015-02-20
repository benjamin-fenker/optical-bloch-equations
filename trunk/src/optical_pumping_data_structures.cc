// Authors: Benjamin Fenker 2013
// Copyright 2012 Benjamin Fenker

#include <cstring>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#include <gsl/gsl_const_mksa.h>
#include "optical_pumping_data_structures.h"
#include "units.h"


extern bool op_verbose;
bool isZero;

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
  bool debug = false;
  double I_s = M_PI / (3.0*tau);
  I_s /= pow(_speed_of_light, 2.0);
  I_s *= _planck_h * pow(nu, 3);
  saturation_intensity = I_s;
  //  saturation_intensity = 1.7397808 * _mW/_cm/_cm;
  if (debug) printf("I_s = %10.8G mW/cm^2\n", saturation_intensity/(_mW/_cm2));
}

void Laser_data::set_field_components() {
  // field[0] = sigma^- and field[1] = sigma^+ are the two basis vectors.
  // No capability at this time to do x-hat and y-hat basis
  if (stokes[0] < 0.0) {
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
  if (field[0] < 0.0 || field[2] < 0.0) {
    printf("Electric field not possible.  Aborting\n");
    exit(1);
  }
  for (int i = 0; i < 3; i++) {
    intensity[i] = 0.5 * _epsilon_0 * _speed_of_light * pow(field[i], 2.0);
    // printf("intensity[%d] = %8.6G mW/cm2\n", i, intensity[i]/(_mW/_cm2));
  }
}

void Laser_data::switch_off(double tau) {
  power = 0.0;
  set_saturation_intensity(tau);
  for (int i = 0; i < 4; i++) {
    stokes[i] = 0.0;
  }
  set_field_components();
  set_intensity_components();
}

DM_container::DM_container(int setNumEStates, int setNumFStates,
                           int setNumGStates) :
    numEStates(setNumEStates), numFStates(setNumFStates),
    numGStates(setNumGStates),
    ee(numEStates,
       vector<gsl_complex>(numEStates, gsl_complex_rect(0.0, 0.0))),
    ff(numFStates,
       vector<gsl_complex>(numFStates, gsl_complex_rect(0.0, 0.0))),
    gg(numGStates,
       vector<gsl_complex>(numGStates, gsl_complex_rect(0.0, 0.0))),
    ef(numEStates,
       vector<gsl_complex>(numFStates, gsl_complex_rect(0.0, 0.0))),
    eg(numEStates,
       vector<gsl_complex>(numGStates, gsl_complex_rect(0.0, 0.0))),
    fg(numFStates,
       vector<gsl_complex>(numGStates, gsl_complex_rect(0.0, 0.0))) {
}

void DM_container::add(DM_container* dm, DM_container *other) {
  for (int e = 0; e < dm -> numEStates; e++) {
    for (int ep = 0; ep < dm -> numEStates; ep++) {
      dm -> ee[e][ep] = gsl_complex_add(dm -> ee[e][ep],
                                            other -> ee[e][ep]);
    }
    for (int f = 0; f < dm -> numFStates; f++) {
      dm -> ef[e][f] = gsl_complex_add(dm -> ef[e][f],
                                           other -> ef[e][f]);
    }
    for (int g = 0; g < dm -> numGStates; g++) {
      dm -> eg[e][g] = gsl_complex_add(dm -> eg[e][g],
                                           other -> eg[e][g]);
    }
  }
  for (int f = 0; f < dm -> numFStates; f++) {
    for (int fp = 0; fp < dm -> numFStates; fp++) {
      dm -> ff[f][fp] = gsl_complex_add(dm -> ff[f][fp],
                                            other -> ff[f][fp]);
    }
    for (int g = 0; g < dm -> numGStates; g++) {
      dm -> fg[f][g] = gsl_complex_add(dm -> fg[f][g],
                                           other -> fg[f][g]);
    }
  }
  for (int g = 0; g < dm -> numGStates; g++) {
    for (int gp = 0; gp < dm -> numGStates; gp++) {
      dm -> gg[g][gp] = gsl_complex_add(dm -> gg[g][gp],
                                            other -> gg[g][gp]);
    }
  }
}

bool DM_container::okay_to_add(DM_container *dm, DM_container *other) {
  DM_container *tmp = new DM_container(*dm);
  DM_container::add(tmp, other);
  for (int e = 0; e < dm -> numEStates; e++) {
    double epop = GSL_REAL(dm->ee[e][e]);
    if (epop > 1 || epop < 0) return false;
  }
  for (int f = 0; f < dm -> numFStates; f++) {
    double fpop = GSL_REAL(dm->ff[f][f]);
    if (fpop > 1 || fpop < 0) return false;
  }
  for (int g = 0; g < dm -> numGStates; g++) {
    double gpop = GSL_REAL(dm->gg[g][g]);
    if (gpop > 1 || gpop < 0) return false;
  }
  return true;
}

void DM_container::mul(DM_container *dm, double c) {
  for (int e = 0; e < dm -> numEStates; e++) {
    for (int ep = 0; ep < dm -> numEStates; ep++) {
      dm -> ee[e][ep] = gsl_complex_mul_real(dm -> ee[e][ep], c);
    }
    for (int f = 0; f < dm -> numFStates; f++) {
      dm -> ef[e][f] = gsl_complex_mul_real(dm -> ef[e][f], c);
    }
    for (int g = 0; g < dm -> numGStates; g++) {
      dm -> eg[e][g] = gsl_complex_mul_real(dm -> eg[e][g], c);
    }
  }
  for (int f = 0; f < dm -> numFStates; f++) {
    for (int fp = 0; fp < dm -> numFStates; fp++) {
      dm -> ff[f][fp] = gsl_complex_mul_real(dm -> ff[f][fp], c);
    }
    for (int g = 0; g < dm -> numGStates; g++) {
      dm -> fg[f][g] = gsl_complex_mul_real(dm -> fg[f][g], c);
    }
  }
  for (int g = 0; g < dm -> numGStates; g++) {
    for (int gp = 0; gp < dm -> numGStates; gp++) {
      dm -> gg[g][gp] = gsl_complex_mul_real(dm -> gg[g][gp], c);
    }
  }
}

bool DM_container::equalsZero() {
  isZero = false;
  double eps = pow(10, -6);
  for (int e = 0; e < numEStates; e++) {
    for (int ep = 0; ep < numEStates; ep++) {
      if (fabs(GSL_REAL(ee[e][ep])) > eps) return false;
      if (fabs(GSL_IMAG(ee[e][ep])) > eps) return false;
    }
    for (int f = 0; f < numFStates; f++) {
      if (fabs(GSL_REAL(ef[e][f])) > eps) return false;
      if (fabs(GSL_IMAG(ef[e][f])) > eps) return false;
    }
    for (int g = 0; g < numGStates; g++) {
      if (fabs(GSL_REAL(eg[e][g])) > eps) return false;
      if (fabs(GSL_IMAG(eg[e][g])) > eps) return false;
    }
  }
  for (int f = 0; f < numFStates; f++) {
    for (int fp = 0; fp < numFStates; fp++) {
      if (fabs(GSL_REAL(ff[f][fp])) > eps) return false;
      if (fabs(GSL_IMAG(ff[f][fp])) > eps) return false;
    }
    for (int g = 0; g < numGStates; g++) {
      if (fabs(GSL_REAL(fg[f][g])) > eps) return false;
      if (fabs(GSL_IMAG(fg[f][g])) > eps) return false;
    }
  }
  for (int g = 0; g < numGStates; g++) {
    for (int gp = 0; gp < numGStates; gp++) {
      if (fabs(GSL_REAL(gg[g][gp])) > eps) return false;
      if (fabs(GSL_IMAG(gg[g][gp])) > eps) return false;
    }
  }
  isZero = true;
  return true;
}

op_parameters::op_parameters() {
  // ****These are only defaults****
  isotope = "37K";
  method = "O";
  out_file = "opData.dat";
  Je2 = 1;
  tune_fe = 4;
  tune_ge = 4;
  tmax = 1.0 * _ns;
  tstep = 1.0 * _ns;
  zeeman = true;
  hyperfine_gr = true;
  hyperfine_ex = true;

  laser_fe.power = 200.0 * (_uW/_cm2);
  laser_fe.detune = -4.5 * _MHz;
  laser_fe.linewidth = 0.2 * _MHz;
  laser_fe.s3s0 = 1.0;
  laser_fe.offtime = -1.0 * _ns;

  laser_ge.power = 200.0 * (_uW/_cm2);
  laser_ge.detune = -4.5 * _MHz;
  laser_ge.linewidth = 0.2 * _MHz;
  laser_ge.s3s0 = 1.0;
  laser_ge.offtime = -1.0 * _ns;

  Bx = 0.0 * _G;
  Bz = 2.0 * _G;

  population_tilt = 0.0;                // obsolete

  verbosity = 0;

  rf_linewidth = 500.0 * _Hz;
  initial_population = "uniform";
}

op_parameters::op_parameters(int argc, char* argv[]) {
  isotope = argv[1];
  if (argc > 2) Je2 = atoi(argv[2]);
  if (argc > 3) tmax = atof(argv[3])*_ns;
  if (argc > 4) tstep = atof(argv[4])*_ns;
  if (argc > 5) Bz = atof(argv[5])*_G;
  if (argc > 6) Bx = atof(argv[6])*_G;
  if (argc > 7) laser_fe.detune = atof(argv[7]) *_MHz;
  if (argc > 8) laser_ge.detune = atof(argv[8]) * _MHz;
  if (argc > 9) {
    laser_fe.linewidth = atof(argv[9]) *_MHz;
    laser_ge.linewidth = atof(argv[9]) *_MHz;
  }
  if (argc > 10) laser_fe.power = atof(argv[10]) * _uW/_cm2;
  if (argc > 11) laser_ge.power = atof(argv[11]) * _uW/_cm2;
  if (argc > 12) laser_fe.s3s0 = atof(argv[12]);
  if (argc > 13) laser_ge.s3s0 = atof(argv[13]);
  if (argc > 14) out_file = string(argv[14]);
  if (argc > 15) method = string(argv[15]);
  if (argc > 16) population_tilt = atof(argv[16]);
}

op_parameters::op_parameters(FILE *file) {
  char expectedInput[40] = "file";
  string file_s = "tmp.dat";
  readAndCheckFromFile(file, expectedInput, &file_s);
  out_file = file_s;
  

  snprintf(expectedInput, sizeof(expectedInput), "method");
  readAndCheckFromFile(file, expectedInput, &method);

  snprintf(expectedInput, sizeof(expectedInput), "isotope");
  readAndCheckFromFile(file, expectedInput, &isotope);

  snprintf(expectedInput, sizeof(expectedInput), "Je2");
  readAndCheckFromFile(file, expectedInput, &Je2);
  
  snprintf(expectedInput, sizeof(expectedInput), "tstep");
  readAndCheckFromFile(file, expectedInput, &tstep);
  tstep *= _ns;              // Have to get the units right!

  snprintf(expectedInput, sizeof(expectedInput), "tmax");
  readAndCheckFromFile(file, expectedInput, &tmax);
  tmax *= _ns;              // Have to get the units right!

  int useCoherence = 1;
  snprintf(expectedInput, sizeof(expectedInput), "zCoherences");
  readAndCheckFromFile(file, expectedInput, &useCoherence);
  if (useCoherence != 1) {
    zeeman = false;
  } else {
    zeeman = true;
  }

  useCoherence = 1;
  snprintf(expectedInput, sizeof(expectedInput), "hfCoherences_gr");
  readAndCheckFromFile(file, expectedInput, &useCoherence);
  if (useCoherence != 1) {
    hyperfine_gr = false;
  } else {
    hyperfine_gr = true;
  }

  useCoherence = 1;
  snprintf(expectedInput, sizeof(expectedInput), "hfCoherences_ex");
  readAndCheckFromFile(file, expectedInput, &useCoherence);
  if (useCoherence != 1) {
    hyperfine_ex = false;
  } else {
    hyperfine_ex = true;
  }

  snprintf(expectedInput, sizeof(expectedInput), "laser_fe_power");
  readAndCheckFromFile(file, expectedInput, &laser_fe.power);
  laser_fe.power *= _uW/_cm2;            // Have to get the units right!

  snprintf(expectedInput, sizeof(expectedInput), "laser_fe_s3");
  readAndCheckFromFile(file, expectedInput, &laser_fe.s3s0);

  snprintf(expectedInput, sizeof(expectedInput), "laser_fe_linewidth");
  readAndCheckFromFile(file, expectedInput, &laser_fe.linewidth);
  laser_fe.linewidth *= _MHz;            // Have to get the units right!

  snprintf(expectedInput, sizeof(expectedInput), "laser_fe_nomTune");
  readAndCheckFromFile(file, expectedInput, &tune_fe);

  snprintf(expectedInput, sizeof(expectedInput), "laser_fe_detune");
  readAndCheckFromFile(file, expectedInput, &laser_fe.detune);
  laser_fe.detune *= _MHz;            // Have to get the units right!

  snprintf(expectedInput, sizeof(expectedInput), "laser_fe_offTime");
  readAndCheckFromFile(file, expectedInput, &laser_fe.offtime);
  laser_fe.offtime *= _ns;            // Have to get the units right!

  snprintf(expectedInput, sizeof(expectedInput), "laser_ge_power");
  readAndCheckFromFile(file, expectedInput, &laser_ge.power);
  laser_ge.power *= _uW/_cm2;            // Have to get the units right!

  snprintf(expectedInput, sizeof(expectedInput), "laser_ge_s3");
  readAndCheckFromFile(file, expectedInput, &laser_ge.s3s0);

  snprintf(expectedInput, sizeof(expectedInput), "laser_ge_linewidth");
  readAndCheckFromFile(file, expectedInput, &laser_ge.linewidth);
  laser_ge.linewidth *= _MHz;            // Have to get the units right!

  snprintf(expectedInput, sizeof(expectedInput), "laser_ge_nomTune");
  readAndCheckFromFile(file, expectedInput, &tune_ge);

  snprintf(expectedInput, sizeof(expectedInput), "laser_ge_detune");
  readAndCheckFromFile(file, expectedInput, &laser_ge.detune);
  laser_ge.detune *= _MHz;            // Have to get the units right!

  snprintf(expectedInput, sizeof(expectedInput), "laser_ge_offTime");
  readAndCheckFromFile(file, expectedInput, &laser_ge.offtime);
  laser_ge.offtime *= _ns;            // Have to get the units right!

  snprintf(expectedInput, sizeof(expectedInput), "B_z");
  readAndCheckFromFile(file, expectedInput, &Bz);
  Bz *= _G;            // Have to get the units right!

  snprintf(expectedInput, sizeof(expectedInput), "B_x");
  readAndCheckFromFile(file, expectedInput, &Bx);
  Bx *= _G;            // Have to get the units right!

  snprintf(expectedInput, sizeof(expectedInput), "tilt");
  readAndCheckFromFile(file, expectedInput, &population_tilt);

  snprintf(expectedInput, sizeof(expectedInput), "rf_linewidth");
  readAndCheckFromFile(file, expectedInput, &rf_linewidth);
  rf_linewidth *= _Hz;

  snprintf(expectedInput, sizeof(expectedInput), "population");
  readAndCheckFromFile(file, expectedInput, &initial_population);
  //  printf("Read %s\n", population_setup.c_str());
}

void op_parameters::readAndCheckFromFile(FILE *f, char *parameter, string *s) {
  char tempL[40] = "";
  char tempR[40] = "";
  int val;
  val = fscanf(f, "%s\t%s", tempL, tempR);
  if (val != 2) {
    printf("Read error.\n");
    exit(1);
  }
  // printf("READ %s\n", tempR);
  if (strcmp(parameter, tempL) != 0) {
    printf("Unexpected line in input file: %s\n",
           tempR);
    exit(1);
  }
  string out = tempR;
  *s = out;
}

void op_parameters::readAndCheckFromFile(FILE *f, char *parameter, int *i) {
  char temp[20] = "";
  int val = fscanf(f, "%s\t%d", temp, i);
  if (val != 2) {
    printf("Read error.\n");
    exit(1);
  }
  if (strcmp(parameter, temp) != 0) {
    printf("Unexpected line in input file: %s\n",
           temp);
    exit(1);
  }
}

void op_parameters::readAndCheckFromFile(FILE *f, char *parameter, double *i) {
  char temp[20] = "";
  int val = fscanf(f, "%s\t%lf", temp, i);
  if (val != 2) {
    printf("Read error.\n");
    exit(1);
  }
  if (strcmp(parameter, temp) != 0) {
    printf("Unexpected line in input file: %s\n",
           temp);
    exit(1);
  }
}

void op_parameters::PrintToFile(string fname) {
  FILE *file;
  file = fopen(fname.c_str(), "w");
  
  fprintf(file, "file %s\n", out_file.c_str());
  fprintf(file, "method %s\n", method.c_str());
  fprintf(file, "isotope %s\n", isotope.c_str());
  fprintf(file, "Je2 %d\n", Je2);
  fprintf(file, "tstep %f\n", tstep/_ns);
  fprintf(file, "tmax %f\n", tmax/_ns);
  fprintf(file, "zCoherences %d\n", (zeeman)? 1:0);
  fprintf(file, "hfCoherences_gr %d\n", (hyperfine_gr)? 1:0);
  fprintf(file, "hfCoherences_ex %d\n", (hyperfine_ex)? 1:0);
  fprintf(file, "laser_fe_power %f\n", laser_fe.power/(_uW/_cm2));
  fprintf(file, "laser_fe_s3 %f\n", laser_fe.s3s0);
  fprintf(file, "laser_fe_linewidth %f\n", laser_fe.linewidth/_MHz);
  fprintf(file, "laser_fe_nomTune %d\n", tune_fe);
  fprintf(file, "laser_fe_detune %f\n", laser_fe.detune/_MHz);
  fprintf(file, "laser_fe_offTime %f\n", laser_fe.offtime/_ns);
  fprintf(file, "laser_ge_power %f\n", laser_ge.power/(_uW/_cm2));
  fprintf(file, "laser_ge_s3 %f\n", laser_ge.s3s0);
  fprintf(file, "laser_ge_linewidth %f\n", laser_ge.linewidth/_MHz);
  fprintf(file, "laser_ge_nomTune %d\n", tune_ge);
  fprintf(file, "laser_ge_detune %f\n", laser_ge.detune/_MHz);
  fprintf(file, "laser_ge_offTime %f\n", laser_ge.offtime/_ns);
  fprintf(file, "B_z %f\n", Bz/_G);
  fprintf(file, "B_x %f\n", Bx/_G);
  fprintf(file, "tilt %f\n", population_tilt);
  fprintf(file, "rf_linewidth %f\n", rf_linewidth/_Hz);
  fprintf(file, "population %s\n", initial_population.c_str());
  fprintf(file, "\n");

  fprintf(file, "Units must be ns, uW/cm^2, MHz, G where appropriate\n\n");
  fprintf(file, "Laser offTime should be < 0 for 'always on' (the usual case)\n\n");
  fprintf(file, "Each 'nomTune' is the 2*Mf of the state\n\n");
  fprintf(file, "It is always assumed that your are pumping to the HIGHEST F-manifold\n");
  fprintf(file, "in the excited state\n\n");
  fprintf(file, "For rf_linewidth the UNITS ARE HERTZ. If rf_linewidth set to >= 0.0,\n");
  fprintf(file, "the (Gam1+Gam2)/2 term in Eq 37 will replaced by this value.  This\n");
  fprintf(file, "simulates 100%% correlated lasers as in Gu (2003).  If rf_linewidth < \n");
  fprintf(file, "0, use Eq 37 (T&J) exactly as is, simulating 0%% correlated lasers. xb\n");

  fclose(file);
}
