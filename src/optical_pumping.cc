// Authors: Benjamin Fenker 2013
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
// extern char outFile[50];
extern bool op_batch;
extern bool isZero;

int OpticalPumping::pump(op_parameters params) {
  Laser_parameters fe = params.laser_fe;
  Laser_parameters ge = params.laser_ge;

  int status = pump(params.isotope, params.method, params.tmax, params.tstep,
                    params.zeeman, params.hyperfine_ex, params.hyperfine_gr,
                    params.Je2, params.tune_fe, params.tune_ge, fe.power, ge.power,
                    fe.detune, ge.detune, fe.linewidth, ge.linewidth, fe.s3s0, ge.s3s0,
                    fe.offtime, ge.offtime, params.Bz, params.Bx, params.out_file,
                    params.population_tilt, params.verbosity, params.rf_linewidth);
  return status;
}
int OpticalPumping::pump(string isotope, string method, double tmax,
                         double tStep, bool zCoherences, bool hfCoherences_ex,
                         bool hfCoherences_gr, int temp_Je2,
                         int nominalSublevelTune2_ef,
                         int nominalSublevelTune2_eg,
                         double laser_fe_power, double laser_ge_power,
                         double laser_fe_detune, double laser_ge_detune,
                         double laser_fe_linewidth, double laser_ge_linewidth,
                         double laser_fe_s3_over_s0, double laser_ge_s3_over_s0,
                         double laser_fe_offTime, double laser_ge_offTime,
                         double set_B_z, double set_B_x, string outFile,
                         double tilt, int verbose, double rf_linewidth) {
  bool debug = false;
  atom_data atom;
  atom.Je2 = temp_Je2;
  isZero = false;
  magnetic_field_data field;
  // **Physical constant**
  // printf("mu_B = %10.8G MHz/G\n", _bohr_magneton/_planck_h/(_MHz/_G));
  field.B_z = set_B_z;
  field.B_x = set_B_x;

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
  atom.linewidth = 1.0 / (2*M_PI*atom.tau);

  // These will hold the F quantum number of the two ground and one excited
  // state that the lasers will be detuned from
  int tuned_F_F2 = atom.I2+1;
  int tuned_G_F2 = abs(atom.I2-1);
  int tuned_E_F2 = atom.Je2 + atom.I2;
  // *************************************************************************

  double laser_fe_nu = alk.getLaserFrequency(atom, field, tuned_F_F2,
                                             nominalSublevelTune2_ef,
                                             tuned_E_F2,
                                             //nominalSublevelTune2_ef,
                                             // tuned_E_F2,
                                             nominalSublevelTune2_ef,
                                             // Tunes max sublevel
                                             laser_fe_detune);

  // double laser_fe_nu = alk.getLaserFrequency(atom, field, 4,
  //                                            4,
  //                                            4,
  //                                            4,
  //                                            // tuned_E_F2,
  //                                            // Tunes max sublevel
  //                                            laser_fe_detune);

  double laser_ge_nu = alk.getLaserFrequency(atom, field, tuned_G_F2,
                                             nominalSublevelTune2_eg,
                                             tuned_E_F2,
                                             nominalSublevelTune2_eg,
                                             // tuned_E_F2,
                                             //0,
                                             // Tunes max sublevel
                                             laser_ge_detune);
  // laser_fe_nu = _speed_of_light * 12985. / _cm;
  // laser_ge_nu = _speed_of_light * 12985. / _cm;

  // printf("Laser fe: %f MHz\n", laser_fe_nu/_MHz);
  // printf("Laser ge: %f MHz\n", laser_ge_nu/_MHz);
  
  Laser_data laser_fe(laser_fe_nu, laser_fe_power, laser_fe_detune,
                      laser_fe_linewidth, laser_fe_s3_over_s0, atom.tau);
  Laser_data laser_ge(laser_ge_nu, laser_ge_power, laser_ge_detune,
                      laser_ge_linewidth, laser_ge_s3_over_s0, atom.tau);

  // Same units as atom.Aj_g,e
  atom.gamma_spon = (alk.getGamma(atom.tau, laser_fe.power));

  string longMethod;
  if (method == "O") {
    longMethod = "Optical-Bloch Equations";
  } else if (method == "R") {
    longMethod = "Rate Equations";
  } else {
    printf("Method not found.  Use \"R\" or \"O\"\n");
    return 2;
  }

  if (verbose >= 0) {
    printf("** OPTICAL PUMPING **\n");
    printf("%s\n\n", longMethod.c_str());

    printf("Pumping %s\n\tAs = %8.6G MHz \n ", isotope.c_str(), atom.Aj_g/_MHz);
    printf("\tAp = %8.6G MHz \n \tg-Factor = %8.6G\n", atom.Aj_e/_MHz,
           atom.g_I);
    printf("\tatom.nu_excited = %14.10G MHz\n", atom.nu_excited/_MHz);
    if (atom.I2%2 == 0) {
      printf("\tNuclear spin: %d \n", atom.I2/2);
    } else {
      printf("\tNuclear spin: %d/2 \n", atom.I2);
    }
    printf("\tG States: %i\tF States: %i\tE States %i\n\n",
           atom.numGStates, atom.numFStates, atom.numEStates);


    printf("\tTau: %5.2G ns\n\tgamma_spon = %5.2G MHz \n\n", atom.tau/_ns,
           atom.gamma_spon/_MHz);
    printf("Tmax: %8.1G ns \nTime Step: %4.2G ns \n", tmax/_ns, tStep/_ns);
    printf("Magnetic Field: %+6.4G z + %+6.4G x G\n\n", field.B_z/_G,
           field.B_x/_G);

    printf("\nLaser g->e (Laser 1) data:\n\tDetuned %+5.2G MHz from the |",
           laser_ge.detune/_MHz);
    printf("%i/2,%i/2> ---> |%i/2,%i/2>", tuned_G_F2, nominalSublevelTune2_eg,
           tuned_E_F2, nominalSublevelTune2_eg);
    printf(" transition.\n\tFrequency = %14.10G MHz\n\tLinewidth = %5.2G MHz\n",
           laser_ge.nu/_MHz, laser_ge.linewidth/_MHz);
    printf("\tIntensity = %5.2G mW/cm^2\n", laser_ge.power/(_mW/_cm2));
    printf("\ts3 = %+5.3G V^2/m^2--> I = <%8.6G, %8.6G> mW/cm^2\n",
           laser_ge.stokes[3]/(_V*_V/_m*_m), laser_ge.intensity[0]/(_mW/_cm2),
           laser_ge.intensity[2]/(_mW/_cm2));
    printf("\t                         E = <%8.6G, %8.6G> V/m\n",
           laser_ge.field[0]/(_V/_m), laser_ge.field[2]/(_V/_m));

    printf("\nLaser f->e (Laser 2) data:\n\tDetuned %+5.2G MHz from the |",
           laser_fe.detune/_MHz);
    printf("%i/2,%i/2> ---> |%i/2,%i/2>", tuned_F_F2, nominalSublevelTune2_ef,
           tuned_E_F2, nominalSublevelTune2_ef);
    printf(" transition.\n\tFrequency = %14.10G MHz\n\tLinewidth = %5.2G MHz\n",
           laser_fe.nu/_MHz, laser_fe.linewidth/_MHz);
    printf("\tIntensity = %5.2G mW/cm^2\n", laser_fe.power/(_mW/_cm2));
    printf("\ts3 = %+5.3G V^2/m^2 --> I = <%8.6G, %8.6G> mW/cm^2\n",
           laser_fe.stokes[3]/(_V*_V/_m*_m), laser_fe.intensity[0]/(_mW/_cm2),
           laser_fe.intensity[2]/(_mW/_cm2));
    printf("\t                          E = <%8.6G, %8.6G> V/m\n",
           laser_fe.field[0]/(_V/_m), laser_fe.field[2]/(_V/_m));
    if (rf_linewidth >= 0) {
      printf("Lasers are 100%% correlated with RF linewidth %g Hz\n\n", rf_linewidth/_Hz);
    } else {
      printf("Lasers are 0%% correlated\n\n");
    }
  } else if (verbose == -1) {
    printf("Isotope: %s     s_3 = %+6.4f/%6.4f     B_x = %6.1f mG\n",
           isotope.c_str(), laser_fe_s3_over_s0, laser_ge_s3_over_s0,
           field.B_x/_mG);
  }

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

  // Maximum number of time steps to write to a file
  const int max_out = 5000;
  double print_frequency;  // Number of ns between print statements
  //  int total_print;  // Total number of print statements
  print_frequency = tStep;
  while(tmax  / print_frequency > max_out) {
    print_frequency *= 2.0;
  }

  if(op_batch) print_frequency = tStep;
  if (verbose >= 0) printf("Print frequency: %g\n", print_frequency);
  // print_frequency = 10*_us;
  // Setup File I/O for later use

  FILE * file;
  file = fopen(outFile.c_str(), "w");
  if (!debug) {
    if (verbose >= 0) printf("Opening file...%s\n", outFile.c_str());
    if (file == NULL) {
      printf("could not open file %s\n", outFile.c_str());
    }
    // This section is commented out to avoid printing the header info
    /*
    if (!op_batch) {
      // Print data for plotting
      fprintf(file, "%d \t ", atom.numEStates+atom.numFStates+atom.numGStates);
      fprintf(file, "%d \t %6.2G \t %d \t %d \t ", atom.numEStates, 1.0,
              total_print, atom.I2);
      fprintf(file, "%d \t 4.0\n", atom.Je2);

      // Print laser data
      fprintf(file, "%4.2G \t %4.2G \t %4.2G \t %4.2G \t %d \t %4.2G \n ",
              laser_fe.power/(_mW/_cm2), laser_fe.nu/_MHz,
              laser_fe.detune/_MHz, laser_fe.stokes[3]/laser_fe.stokes[0],
              nominalSublevelTune2_ef, laser_fe_offTime/_us);
      fprintf(file, "%4.2G \t %4.2G \t %4.2G \t %4.2G \t %d \t %4.2G \n ",
              laser_ge.power/(_mW/_cm2), laser_ge.nu/_MHz,
              laser_ge.detune/_MHz, laser_ge.stokes[3]/laser_ge.stokes[0],
              nominalSublevelTune2_eg, laser_ge_offTime/_us);

      // Print magnetic field data
      fprintf(file, "%4.2G \t %4.2G \n", field.B_z/_G, field.B_x/_G);

      // Print atomic data
      if (zCoherences) {
        fprintf(file, "%d \t ", 1);
      } else {
        fprintf(file, "%d \t ", 0);
      }
      if (hfCoherences_gr) {
        fprintf(file, "%d \t ", 1);
      } else {
        fprintf(file, "%d \t ", 0);
      }
      if (hfCoherences_ex) {
        fprintf(file, "%d \n ", 1);
      } else {
        fprintf(file, "%d \n ", 0);
      }
      fprintf(file, "%s\n", isotope.c_str());
      fprintf(file, "%s\n", method.c_str());
    }
    */
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
    equ = new Rate_Equations(decomp, laser_fe, laser_ge, tilt);
    //    update_func = &Rate_Equations::update_population_gsl;
    //    int (*update_func)(double, const double[], double[], void*) =
    //      &Rate_Equations::update_population_gsl;
  } else if (strcmp(method.c_str(), "O") == 0) {
    equ = new Density_Matrix(decomp, laser_fe, laser_ge, flags, tilt);
    //  update_func = &Density_Matrix::update_population_gsl;
  } else {
    printf("Failed to initialize the optical pumping method\t");
    printf("Method = %s.  Aborting.", method.c_str());
    return(2);
  }
  equ -> rf_linewidth = rf_linewidth;

  // double tStep_in = tStep;
  // if (op_batch) tStep = 0.01*_us;
  double time = 0.0;
  double nextPrint = 0.0;
  double nextUpdate = 0.0;
  // printf("updateFreq = %8.6G\n", updateFreq);
  bool laser_fe_Off = false;
  bool laser_ge_Off = false;
  bool print_zero = false;
  int print_progress_width = 60;
  double print_progress_dt = tmax / print_progress_width;
  double next_print_progress = 0.0;
  if (verbose >= -1) {
    for (int i = 0; i < print_progress_width; i++) printf("=");
    printf("\n");
  }
  while (time < tmax) {
    if (!op_batch) {
      if ((fabs(time - nextUpdate))/_ns < pow(10, -2)) {
        //        printf(" t = %8.6G ns\n", time/_ns);
         nextUpdate += updateFreq;
       }
     }
     //    printf("%g\n", time/_ns);
     if ((fabs(time - next_print_progress)) < tStep/2. && !print_zero) {
       //       printf("A\n");

       next_print_progress += print_progress_dt;
       if (verbose >= -1) {
         printf("=");
         fflush(stdout);
       }
    }
    
    // if ((fabs(time - nextPrint))/_ns < pow(10, -2) {
    if ((fabs(time - nextPrint)/_ns < pow(10, -2))) {
      // equ -> print_data(stdout, time);
      equ->print_data(file, time);
      //      equ->print_density_matrix(stdout);
      nextPrint += print_frequency;
    }
    if (!laser_ge_Off && laser_ge_offTime > 0 && time >= laser_ge_offTime) {
      laser_ge_Off = true;
      printf("Switching off ge laser at time = %6.4G ns\n", time/_ns);
      equ -> switch_off_laser(1);
    }
    if (!laser_fe_Off && laser_fe_offTime > 0 && time >= laser_fe_offTime) {
      laser_fe_Off = true;
      printf("Switching off fe laser at time = %6.4G ns\n", time/_ns);
      equ -> switch_off_laser(2);
    }
    if (time > -2.0*_ns) {
      // equ -> change_magnetic_field(0.0);
    }

    if (!isZero) {
      // **********************************************************************
      equ->update_population_RK4(tStep);  // ********************************
      //      equ -> update_population_euler(tStep);
      // **********************************************************************
    } else {
      if (!print_zero) {
        print_zero = true;
        if (verbose >= -1) {
          printf("Reached steady-state at t = %g us\n", time/_us);
        }
      }
    }
    //   equ->print_density_matrix(stdout);
    if (!equ->is_hermitian()) {
      printf("DENSITY MATRIX NOT HERMITIAN AT t = %4.2G ns\n", time/_ns);
      if (verbose >= 0) {
        equ -> print_data(stdout, time);
        equ -> print_data(file, time);
      }
      return 1;
    }
    if (fabs(equ->get_total_population() - 1.0) > pow(10, -3)) {
      printf("PARTICLES NOT CONSERVED AT t = %4.2G ns\n", time/_ns);
      if (verbose >= 0) {
        equ -> print_data(stdout, time);
        equ -> print_data(file, time);      
      }
      return 1;
    }
    time += tStep;
    // if (time >= 3.999*_us) {
    //   tStep = tStep_in;
    //   print_frequency = 0.2*_us;
    // }
  }
  fclose(file);
  if (verbose >= -1) printf("\n");
  if (verbose >= 0) printf("Output complete and closed: %s\n", outFile.c_str());
  //  equ -> print_density_matrix(stdout);
  delete equ;
  return 0;
}

int OpticalPumping::test() {
  printf("here rwar\n");
  return -1;
}
