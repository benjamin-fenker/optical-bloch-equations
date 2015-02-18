// Authors: Benjamin Fenker 2014


#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <iomanip>
#include "optical_pumping.h"
#include "optical_pumping_data_structures.h"
#include "units.h"

using std::string;

bool op_verbose = false;
bool op_batch = false;

char outFile[50] = "opData.dat";


int main(int argc, char* argv[]) {
  if (!op_batch) printf("\n");
  if (op_batch && op_verbose) {
    printf("Can't be both batch and verbose. See setup_optical_pumping.cc\n");
    exit(1);
  }

  // *******************************
  op_parameters *p_params = new op_parameters();
  // Also will accept command line input
  if (argc > 1) {
    if (strcmp(argv[1], "-h") == 0) {
      printf("Run this program with from an input file or a command line:\n");
      printf("To run from a file do: opticalPumping -f fileName\n\n");
      printf("To run from the command line, give a list of parameters:\n");
      printf("First paramter is isotope (37K, 38K, 41K, 47K) [%s]\n",
             p_params -> isotope.c_str());
      printf("Second parameter is Je2 (1 for D1 line, 3 for D2 line [%d]\n",
             p_params -> Je2);
      printf("Third parameter is tmax in ns [%3.1G]\n", (p_params -> tmax)/_ns);
      printf("Fourth paramter is dt in ns [%3.1G]\n", (p_params -> tstep)/_ns);
      printf("Fifth parameter is B_z in G [%3.1G]\n", (p_params -> Bz)/_G);
      printf("Sixth parameter is B_x in G [%3.1G]\n", (p_params -> Bx)/_G);
      printf("Seventh parameter is f --> e laser detune in MHz [%3.1G]\n",
             (p_params -> laser_fe.detune)/_MHz);
      printf("Eigth parameter is g --> e laser detune in MHz [%3.1G]\n",
             (p_params -> laser_ge.detune)/_MHz);
      printf("Ninth parameter is linewidth (both lasers) in MHz [%3.1G]\n",
             (p_params -> laser_fe.linewidth)/_MHz);
      printf("Tenth parameter is f --> e laser power in uW/cm^2 [%5.3G]\n",
             (p_params -> laser_fe.power)/(_uW/_cm2));
      printf("Eleventh parametr is g --> e laser power in uW/cm^2 [%5.3G]\n",
             (p_params -> laser_ge.power)/(_uW/_cm2));
      printf("Twelth parameter is g --> laser S3 parameter [%3.1G]\n",
             (p_params -> laser_ge.s3s0));
      printf("Thirteenth parameter is f --> laser S3 parameter [%3.1G]\n",
             (p_params -> laser_fe.s3s0));
      printf("Fourteenth parameter is output file [%s]\n",
             p_params -> out_file.c_str());
      printf("Fifteenth parameter is method [%s]\n",
             p_params -> method.c_str());
      printf("Sixteenth parameter is tilt [%g]\n",
             p_params -> population_tilt);
      printf("\n\n");
      return 0;
    } else if (strcmp(argv[1], "-f") == 0) {  // accept input from file
      if (argc == 2) {                        // No file name given
        printf("File name required with -f option.\n");
        printf("opticalPumping -f in.in\n");
        exit(1);
      }
      FILE *file;
      char* inFile = argv[2];
      // printf("Opening file %s\n", argv[2]);
      file = fopen(inFile, "r");
      if (file != NULL) {
        delete p_params;
        p_params = new op_parameters(file);

      } else {                  // File does not exist
        printf("File %s does not exist\n", inFile);
        exit(1);
      }
    } else {
      delete p_params;
      p_params = new op_parameters(argc, argv);
    }
  }

  if (p_params == NULL) p_params = new op_parameters();
  op_parameters params = *p_params;
  OpticalPumping pumper;
  int status = pumper.pump(params);

  if (!op_batch) printf("\nCompleted with status = %d\n\n", status);
  return status;
}
