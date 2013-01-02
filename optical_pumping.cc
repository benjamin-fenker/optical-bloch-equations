// Copyright 2012 Benjamin Fenker

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_errno.h>      // GSL Reserves error codes -2 to 32 (succes = 0)
#include <gsl/gsl_matrix.h>

#include "include/optical_pumping.h"
#include "include/alkali.h"
#include "include/eigenvector_helper.h"
#include "include/optical_pumping_method.h"
#include "include/rate_equations.h"
#include "include/density_matrix.h"
#include "include/optical_pumping_data_structures.h"
#include "include/units.h"

using std::string;
using std::max;

int OpticalPumping::pump(string isotope, string method, double tmax,
                         double tStep, bool zCoherences, bool hfCoherences_ex,
                         bool hfCoherences_gr, int temp_Je2,
                         double laser_fe_power, double laser_ge_power,
                         double laser_fe_detune, double laser_ge_detune,
                         double laser_fe_linewidth, double laser_ge_linewidth,
                         double laser_fe_s3_over_s0, double laser_ge_s3_over_s0,
                         double set_B_z) {
  bool debug = false;
  atom_data atom;
  atom.Je2 = temp_Je2;

  magnetic_field_data field;
  // **Physical constant**
  field.mu_B = 1.399625 * (_MHz/_G);  // MHz/Gauss
  // printf("mu_B = %10.8G MHz/G\n", field.mu_B/(_MHz/_G));
  field.B_z = set_B_z;
  field.B_x = 0.0;

  coherence_flags flags;
  flags.zCoherences = zCoherences;
  flags.hfCoherences_ex = hfCoherences_ex;
  flags.hfCoherences_gr = hfCoherences_gr;

  // Create helper classes to do useful things
  Alkali alk;

  if (alk.lookupParameters(isotope, atom.Je2, &atom.I2, &atom.Aj_g,
                             &atom.Aj_e, &atom.g_I, &atom.nu_excited,
                             &atom.tau)) {
    printf("Parameter lookup failed.\n");
    printf("Isotope = %s \t Je2 = %d\n", isotope.c_str(), temp_Je2);
    return 1;
  }
  atom.numEStates = alk.getNumberOfExcitedStates(atom.I2, atom.Je2);
  atom.numFStates = alk.getNumberOfGroundStates_f(atom.I2);
  atom.numGStates = alk.getNumberOfGroundStates_g(atom.I2);
  atom.numBasisStates_ground = atom.numGStates + atom.numFStates;
  atom.numBasisStates_excited = atom.numEStates;
  atom.linewidth = 1 / atom.tau;

  // These will hold the F quantum number of the two ground and one excited
  // state that the lasers will be detuned from
  int tuned_F_F2 = atom.I2+1;
  int tuned_G_F2 = abs(atom.I2-1);
  int tuned_E_F2 = atom.Je2 + atom.I2;
  // *************************************************************************

  double laser_fe_nu = alk.getLaserFrequency(atom, tuned_F_F2, tuned_E_F2,
                                      laser_fe_detune);
  // Same units as atom.Aj_g,e
  double laser_ge_nu = alk.getLaserFrequency(atom, tuned_G_F2, tuned_E_F2,
                                      laser_ge_detune);
  Laser_data laser_fe(laser_fe_nu, laser_fe_power, laser_fe_detune,
                      laser_fe_linewidth, laser_fe_s3_over_s0, atom.tau);
  Laser_data laser_ge(laser_ge_nu, laser_ge_power, laser_ge_detune,
                      laser_ge_linewidth, laser_ge_s3_over_s0, atom.tau);

  // Same units as atom.Aj_g,e
  atom.gamma_spon = (alk.getGamma(atom.tau, laser_fe.power));

  printf("** OPTICAL PUMPING **\n\n");
  printf("Time step = %4.2G\n", tStep/_ns);
  printf("Pumping %s\t As = %8.6G MHz \t ", isotope.c_str(), atom.Aj_g/_MHz);
  printf("Ap = %8.6G MHz \t g-Factor = %8.6G\n", atom.Aj_e/_MHz, atom.g_I);
  printf("atom.nu_excited = %14.10G MHz\n", atom.nu_excited/_MHz);
  if (atom.I2%2 == 0) {
    printf("Nuclear spin: %d \n", atom.I2/2);
  } else {
    printf("Nuclear spin: %d/2 \n", atom.I2);
  }
  printf("Tau: %5.2G ns\t gamma_spon = %5.2G MHz \n", atom.tau/_ns,
         atom.gamma_spon/_MHz);
  printf("Tmax: %8.1G ns \t Time Step: %4.2G ns \n", tmax/_ns, tStep/_ns);
  printf("Magnetic Field: %4.2G G\n", field.B_z/_G);
  printf("\nLaser g->e data.  Laser 1: Detuned %5.2G MHz from the |",
         laser_ge.detune/_MHz);
  printf("%i/2,%i> ---> |%i/2,%i>", tuned_G_F2, 0, tuned_E_F2, 0);
  printf(" transition.  \nFrequency = %14.10G MHz, Linewidth = %5.2G MHz",
         laser_ge.nu/_MHz, laser_ge.linewidth/_MHz);
  printf(", Intensity = %5.2G mW/cm^2\t", laser_ge.power/(_mW/_cm2));
  printf(", s3 = %5.3G --> I = <%8.6G, %8.6G> mW/cm^2\t", laser_ge.stokes[3],
         laser_ge.intensity[0]/(_mW/_cm2), laser_ge.intensity[1]/(_mW/_cm2));
  printf(", E = <%8.6G, %8.6G, %8.6G> V/m\n", laser_ge.field[0],
         laser_ge.field[1], laser_ge.field[2]);
  //  printf("Field = %8.6G V/m\n", laser_ge.field/(_V/_m));
  printf("\nLaser f->e data.  Laser 1: Detuned %5.2G MHz from the |",
         laser_fe.detune/_MHz);
  printf("%i/2,%i> ---> |%i/2,%i>", tuned_F_F2, 0, tuned_E_F2, 0);
  printf(" transition.  \nFrequency = %14.10G MHz, Linewidth = %5.2G MHz",
         laser_fe.nu/_MHz, laser_fe.linewidth/_MHz);
  printf(", Intensity = %5.2G mW/cm^2\t", laser_fe.power/(_mW/_cm2));
  printf(", s3 = %5.3G --> I = <%8.6G, %8.6G> mW/cm^2\t", laser_fe.stokes[3],
         laser_fe.intensity[0]/(_mW/_cm2), laser_fe.intensity[1]/(_mW/_cm2));
  printf(", E = <%8.6G, %8.6G, %8.6G> V/m\n", laser_fe.field[0],
         laser_fe.field[1], laser_fe.field[2]);
  //  printf("Field = %8.6G V/m\n", laser_fe.field/(_V/_m));
  printf("G States: %i\tF States: %i\tE States %i\n\n",
         atom.numGStates, atom.numFStates, atom.numEStates);

  // Call function to retrieve vectors holding the Iz/Jz decomposotion of the
  // F, Mf eigenstates.
  // const int numEStates = atom.numEStates;
  // const int numGroundStates = atom.numGStates + atom.numFStates;

  /*
  double *ground_IzJz_Decomp = new double[numGroundStates*numGroundStates];
  double *excited_IzJz_Decomp = new double[numEStates*numEStates];
  decomp->diagH(atom, field, 0, 1, atom.Aj_g, &ground_IzJz_Decomp[0]);
  decomp->diagH(atom, field, 1, atom.Je2, atom.Aj_e, &excited_IzJz_Decomp[0]);
  // Test the decomps:
  printf("TESTING DECOMPS:\n");
  for (int i = 0; i < numEStates*numEStates; i++) {
    printf("%8.6G   ", excited_IzJz_Decomp[i]);
    if ((i+1)%numEStates == 0) printf("\n");
  }
  */
  // ***********************************************************************

  // Setup print statements to only print a reasonable number of times
  // ******CAUTION******
  // For very long times, the program will not correctly print out the small
  // amplitude, high frequency oscillations.  This is because the frequency of
  // the oscillation is higher than the frequency of print statements!
  // One thing to do for testing is to increase max out to some huge number
  // allowing the print frquency to be higher, matching the physics frequency
  // ******CAUTION******

  const int max_out = 1000;  // Maximum number of time steps to write to a file
  double print_frequency;  // Number of ns between print statements
  int total_print;  // Total number of print statements
  if (tmax/tStep <= max_out) {
    print_frequency = tStep;
    total_print = static_cast<int>(tmax/tStep);
  } else {
    print_frequency = tmax / static_cast<double>(max_out);
    total_print = max_out;
  }
  // Setup File I/O for later use
  FILE * file;
  if (!debug) {
    file = fopen("opData.dat", "w");
    fprintf(file, "%d \t ", atom.numEStates+atom.numFStates+atom.numGStates);
    fprintf(file, "%d \t %6.2G \t %d \t %d \t ", atom.numEStates, 1.0,
            total_print, atom.I2);
    fprintf(file, "%d \t 4.0\n", atom.Je2);
    fprintf(file, "%s\n", method.c_str());
  } else {
    file = stdout;
  }
  double updateFreq = max(tStep * 1000.0, 1000.0)*_ns;
  updateFreq = 50*_us;
  // Sets frequency that programs informs user that its still working

  /*
  FILE * eeFile = fopen("eeData.dat", "w");
  FILE * ffFile = fopen("ffData.dat", "w");
  FILE * ggFile = fopen("ggData.dat", "w");
  FILE * efFile = fopen("efData.dat", "w");
  FILE * egFile = fopen("egData.dat", "w");
  FILE * fgFile = fopen("fgData.dat", "w");
  */
  // End set up output

  OpticalPumping_Method *equ;
  Eigenvector_Helper decomp(atom, field);
  //  int (*update_func)(double, const double[], double[], void*);

  if (strcmp(method.c_str(), "R") == 0) {
    equ = new Rate_Equations(decomp, laser_fe, laser_ge);
    //    update_func = &Rate_Equations::update_population_gsl;
    //    int (*update_func)(double, const double[], double[], void*) =
    //      &Rate_Equations::update_population_gsl;
  } else if (strcmp(method.c_str(), "O") == 0) {
    equ = new Density_Matrix(decomp, laser_fe, laser_ge, flags);
    //  update_func = &Density_Matrix::update_population_gsl;
  } else {
    printf("Failed to initialize the optical pumping method\t");
    printf("Method = %s.  Aborting.", method.c_str());
    return(2);
  }

  double time = 0.0;
  double nextPrint = 0.0;
  double nextUpdate = 0.0;
  // printf("updateFreq = %8.6G\n", updateFreq);
  while (time < tmax) {
    if ((fabs(time - nextUpdate))/_ns < pow(10, -2)) {
      printf(" t = %8.6G ns\r", time/_ns);
      nextUpdate += updateFreq;
    }
    if ((fabs(time - nextPrint))/_ns < pow(10, -2)) {
      equ->print_data(file, time);
      //  equ->print_density_matrix(stdout);
      nextPrint += print_frequency;
    }
    equ->update_population(tStep);
    if (!equ->is_hermitian()) {
      printf("DENSITY MATRIX NOT HERMITIAN AT t = %4.2G ns\n", time/_ns);
      equ->print_density_matrix(stdout);
      return 1;
    }
    if (fabs(equ->get_total_population() - 1.0) > pow(10, -6)) {
      printf("PARTICLES NOT CONSERVED AT t = %4.2G ns\n", time/_ns);
      equ->print_data(stdout, time);
      return 1;
    }
    time += tStep;
  }
  /*
  printf("\n");
  gsl_odeiv2_system sys = {update_func, NULL, equ->totalTerms, &(equ->data)};
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys,
                                                       gsl_odeiv2_step_rk8pd,
                                                       1e-6, 1e-6, 0);
  double i;
  double t = 0.0;

  fprintf(file, "%8.6G   ", t/_us);
  if (strcmp(method.c_str(), "R") == 0) {
    for (int j = 0; j < equ->totalTerms; j++) {
      fprintf(file, "%5.3G   ", equ->population[j]);
    }
  }
  if (strcmp(method.c_str(), "O") == 0) {
    int index;
    for (int g = 0; g < atom.numGStates; g++) {
      index = g*((atom.numGStates*2)+2);
      fprintf(file, "%5.3G    ", equ->population[index]);
    }
    int terms = atom.numGStates * atom.numGStates * 2;
    for (int f = 0; f < atom.numFStates; f++) {
      index = terms + f*((atom.numFStates*2)+2);
      fprintf(file, "%5.3G    ", equ->population[index]);
    }
    terms += (atom.numFStates * atom.numFStates * 2);
    for (int e = 0; e < atom.numEStates; e++) {
      index = terms + e*((atom.numEStates*2)+2);
      fprintf(file, "%5.3G   ", equ->population[index]);
    }
  }
  fprintf(file, "\n");
  printf("Starting loop\n tStep = %5.3G   tmax = %5.3G\n", tStep, tmax);
  for (i = tStep; i < tmax; i += print_frequency) {
    double ti = i;
    int status = gsl_odeiv2_driver_apply(d, &t, ti, equ->population);
    if (status != GSL_SUCCESS) {
      printf("error, return value = %d\n", status);
      break;
    }
    fprintf(file, "%8.6G   ", t/_us);
    if (strcmp(method.c_str(), "R") == 0) {
      for (int j = 0; j < equ->totalTerms; j++) {
        fprintf(file, "%5.3G   ", equ->population[j]);
      }
    }
    if (strcmp(method.c_str(), "O") == 0) {
      int index;
      for (int g = 0; g < atom.numGStates; g++) {
        index = g*((atom.numGStates*2)+2);
        fprintf(file, "%5.3G    ", equ->population[index]);
      }
      int terms = atom.numGStates * atom.numGStates * 2;
      for (int f = 0; f < atom.numFStates; f++) {
        index = terms + f*((atom.numFStates*2)+2);
        fprintf(file, "%5.3G    ", equ->population[index]);
      }
      terms += (atom.numFStates * atom.numFStates * 2);
      for (int e = 0; e < atom.numEStates; e++) {
        index = terms + e*((atom.numEStates*2)+2);
        fprintf(file, "%5.3G   ", equ->population[index]);
      }
    }
    fprintf(file, "\n");
  }
  gsl_odeiv2_driver_free(d);
  */
  //  delete[] ground_IzJz_Decomp;
  // delete[] excited_IzJz_Decomp;
  return 0;
}


