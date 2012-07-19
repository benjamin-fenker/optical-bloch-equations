// Copyright 2012 Benjamin Fenker

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

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
                         double laser_fe_pol[], double laser_ge_pol[],
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
  Eigenvector_Helper decomp;

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
                      laser_fe_linewidth, laser_fe_pol, atom.tau);
  Laser_data laser_ge(laser_ge_nu, laser_ge_power, laser_ge_detune,
                      laser_ge_linewidth, laser_ge_pol, atom.tau);

  // Same units as atom.Aj_g,e
  atom.gamma_spon = (alk.getGamma(atom.tau, laser_fe.power));

  printf("** OPTICAL PUMPING **\n\n");
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
  printf("Field = %8.6G V/m\n", laser_ge.field/(_V/_m));
  printf("\nLaser f->e data.  Laser 1: Detuned %5.2G MHz from the |",
         laser_fe.detune/_MHz);
  printf("%i/2,%i> ---> |%i/2,%i>", tuned_F_F2, 0, tuned_E_F2, 0);
  printf(" transition.  \nFrequency = %14.10G MHz, Linewidth = %5.2G MHz",
         laser_fe.nu/_MHz, laser_fe.linewidth/_MHz);
  printf(", Intensity = %5.2G mW/cm^2\t", laser_fe.power/(_mW/_cm2));
  printf("Field = %8.6G V/m\n", laser_fe.field/(_V/_m));
  printf("G States: %i\tF States: %i\tE States %i\n\n",
         atom.numGStates, atom.numFStates, atom.numEStates);

  // Call function to retrieve vectors holding the Iz/Jz decomposotion of the
  // F, Mf eigenstates.
  const int numEStates = atom.numEStates;
  const int numGroundStates = atom.numGStates + atom.numFStates;

  double *ground_IzJz_Decomp = new double[numGroundStates*numGroundStates];
  double *excited_IzJz_Decomp = new double[numEStates*numEStates];
  decomp.diagH(atom, field, 0, 1, atom.Aj_g, &ground_IzJz_Decomp[0]);
  decomp.diagH(atom, field, 1, atom.Je2, atom.Aj_e, &excited_IzJz_Decomp[0]);
  // ***********************************************************************

  // Setup print statements to only print a reasonable number of times
  const int max_out = 1000;  // Maximum number of time steps to write to a file
  double print_frequency;  // Number of ns between print statements
  int total_print;  // Total number of print statements
  if (tmax/tStep <= max_out) {
    print_frequency = tStep;
    total_print = static_cast<int>(tmax/tStep);
    printf("(if) total_print = %d\n", total_print);
  } else {
    print_frequency = tmax / static_cast<double>(max_out);
    printf("print_frequency = %6.4G ns\n", print_frequency/_ns);
    total_print = max_out;
    printf("(else)total_print = %d\n", total_print);
  }
  // Setup File I/O for later use
  FILE * file;
  if (!debug) {
    file = fopen("opData.dat", "w");
    fprintf(file, "%d \t ", atom.numEStates+atom.numFStates+atom.numGStates);
    fprintf(file, "%d \t %6.2G \t %d \t %d \t ", atom.numEStates, 1.0,
            total_print, atom.I2);
    fprintf(file, "%d 4.0\n", atom.Je2);
  } else {
    file = stdout;
  }
  double updateFreq = max(tStep * 1000.0, 10000.0);
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

  printf("fe_freq = %14.10G\t ge_freq = %14.10G\n", laser_fe.nu/_MHz,
         laser_ge.nu/_MHz);
  OpticalPumping_Method *equ;
  if (strcmp(method.c_str(), "R") == 0) {
    equ = new Rate_Equations(atom, field, laser_fe, laser_ge);
  } else if (strcmp(method.c_str(), "O") == 0) {
    equ = new Density_Matrix(atom, field, laser_fe, laser_ge, flags);
  } else {
    printf("Failed to initialize the optical pumping method\t");
    printf("Method = %s.  Aborting.", method.c_str());
    return(2);
  }

  double time = 0.0;
  double nextPrint = 0.0;
  double nextUpdate = 0.0;
  while (time < tmax) {
    if ((fabs(time - nextUpdate))/_ns < pow(10, -4)) {
      printf(" t = %8.6G ns\t", time/_ns);
      nextUpdate += updateFreq;
    }
    if ((fabs(time - nextPrint))/_ns < pow(10, -4)) {
      equ->print_data(file, time);
      // equ->print_density_matrix(stdout);
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

  // Testing GSL ODE stuff
  // See first page in GSL manual - section 26.6
  double mu = 10.0;
  int (*update_func)(double, const double[], double[], void*) =
    &Rate_Equations::update_population_gsl;
  //  int (*jacob_func)(double, const double[], double*, double*,
  //                                void*) =
  //  &Rate_Equations::jacobian;
  gsl_odeiv2_system sys = {update_func, NULL, 2, &mu};
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys,
                                                       gsl_odeiv2_step_rk8pd,
                                                       1e-6, 1e-6, 0.0);
  int i;
  double t = 0.0, t1 = 100.0;
  double y[2] = {1.0, 0.0};
  for (i = 1; i < 100; i++) {
    double ti = i*t1/100.0;
    int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
    if (status != GSL_SUCCESS) {
      printf("error, return value = %d\n", status);
      break;
    }
    // printf("%.5e %.5e %.5e\n", t, y[0], y[1]);
  }
  gsl_odeiv2_driver_free(d);
  // Done testing gsl ODE stuff



  delete[] ground_IzJz_Decomp;
  delete[] excited_IzJz_Decomp;
  return 0;
}


