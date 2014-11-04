// Authors: Benjamin Fenker 2014


#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <iomanip>
#include "include/optical_pumping.h"
#include "include/units.h"

using std::string;

bool op_verbose = false;
bool op_batch = false;

char outFile[50] = "opData.dat";

void readAndCheckFromFile(FILE *f, char *parameter, string *s) {
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

void readAndCheckFromFile(FILE *f, char *parameter, int *i) {
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

void readAndCheckFromFile(FILE *f, char *parameter, double *i) {
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

int main(int argc, char* argv[]) {
  if (!op_batch) printf("\n");
  if (op_batch && op_verbose) {
    printf("Can't be both batch and verbose. See setup_optical_pumping.cc\n");
    exit(1);
  }
  string method = "O";  // O for OBE and R for Rate Equations

  // These have to be `true' for the -f option to work correctly
  bool zCoherences = true;
  bool hfCoherences_ex = true;
  bool hfCoherences_gr = true;

  double laser_fe_I = 200.0 * (_uW/_cm2);  // uW/cm^2
  double laser_ge_I = 145.0 * (_uW/_cm2);  // uW/cm^2

  // (I think) that John is reportings I_+ + I_- / I_+ - I_-, which is
  // equivalent to s3/s0 in the notation of Jackson 3rd edition
  double laser_fe_s3_over_s0 = 1.0;
  double laser_ge_s3_over_s0 = 1.0;

  int nominalSublevelTune2_fe = 4;
  int nominalSublevelTune2_ge = 4;

  // ****These are only defaults****
  double tmax = 1.0 * _ns;  // ns
  double dt = 1.0 * _ns;  // ns

  string isotope = "37K";
  int Je2 = 1;  // D1(1) or D2(3) line

  double laser_fe_detune = -4.5 * _MHz;  // MHz
  double laser_ge_detune = -4.5 * _MHz;  // MHz

  double laser_fe_linewidth = 0.2 *_MHz;  // MHz (FWHM)
  double laser_ge_linewidth = 0.2 *_MHz;  // MHz (FWHM)

  double laser_fe_offTime = -125000.0*_ns;  // <0 implies always on
  double laser_ge_offTime = -1;  // <0 implies aways on
  double B_z = 2.0 * _G;  // G
  double B_x = 0.0 * _G;  // G

  double tilt = 0.0;                // To start with some asymmetry
  // *******************************

  // Also will accept command line input
  if (argc > 1) {
    if (strcmp(argv[1], "-h") == 0) {
      printf("Run this program with from an input file or a command line:\n");
      printf("To run from a file do: opticalPumping -f fileName\n\n");
      printf("To run from the command line, give a list of parameters:\n");
      printf("First paramter is isotope (37K, 38K, 41K, 47K) [%s]\n",
             isotope.c_str());
      printf("Second parameter is Je2 (1 for D1 line, 3 for D2 line [%d]\n",
             Je2);
      printf("Third parameter is tmax in ns [%3.1G]\n", tmax/_ns);
      printf("Fourth paramter is dt in ns [%3.1G]\n", dt/_ns);
      printf("Fifth parameter is B_z in G [%3.1G]\n", B_z/_G);
      printf("Sixth parameter is B_x in G [%3.1G]\n", B_x/_G);
      printf("Seventh parameter is f --> e laser detune in MHz [%3.1G]\n",
             laser_fe_detune/_MHz);
      printf("Eigth parameter is g --> e laser detune in MHz [%3.1G]\n",
             laser_ge_detune/_MHz);
      printf("Ninth parameter is linewidth (both lasers) in MHz [%3.1G]\n",
             laser_fe_linewidth/_MHz);
      printf("Tenth parameter is f --> e laser power in uW/cm^2 [%5.3G]\n",
             laser_fe_I/(_uW/_cm2));
      printf("Eleventh parametr is g --> e laser power in uW/cm^2 [%5.3G]\n",
             laser_ge_I/(_uW/_cm2));
      printf("Twelth parameter is g --> laser S3 parameter [%3.1G]\n",
             laser_ge_s3_over_s0);
      printf("Thirteenth parameter is f --> laser S3 parameter [%3.1G]\n",
             laser_fe_s3_over_s0);
      printf("Fourteenth parameter is output file [%s]\n", outFile);
      printf("Fifteenth parameter is method [%s]\n", method.c_str());
      printf("Sixteenth parameter is tilt [%g]\n", tilt);
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
        char expectedInput[40] = "file";
        string file_s = "tmp.dat";
        readAndCheckFromFile(file, expectedInput, &file_s);
        // strcpy(const_cast<char *>(file_s.c_str()), outFile);
        snprintf(outFile, sizeof(outFile), "%s", file_s.c_str());

        snprintf(expectedInput, sizeof(expectedInput), "method");
        readAndCheckFromFile(file, expectedInput, &method);

        snprintf(expectedInput, sizeof(expectedInput), "isotope");
        readAndCheckFromFile(file, expectedInput, &isotope);

        snprintf(expectedInput, sizeof(expectedInput), "Je2");
        readAndCheckFromFile(file, expectedInput, &Je2);

        snprintf(expectedInput, sizeof(expectedInput), "tstep");
        readAndCheckFromFile(file, expectedInput, &dt);
        dt *= _ns;              // Have to get the units right!

        snprintf(expectedInput, sizeof(expectedInput), "tmax");
        readAndCheckFromFile(file, expectedInput, &tmax);
        tmax *= _ns;              // Have to get the units right!

        int useCoherence = 1;
        snprintf(expectedInput, sizeof(expectedInput), "zCoherences");
        readAndCheckFromFile(file, expectedInput, &useCoherence);
        if (useCoherence != 1) zCoherences = false;

        useCoherence = 1;
        snprintf(expectedInput, sizeof(expectedInput), "hfCoherences_gr");
        readAndCheckFromFile(file, expectedInput, &useCoherence);
        if (useCoherence != 1) hfCoherences_gr = false;

        useCoherence = 1;
        snprintf(expectedInput, sizeof(expectedInput), "hfCoherences_ex");
        readAndCheckFromFile(file, expectedInput, &useCoherence);
        if (useCoherence != 1) hfCoherences_ex = false;

        snprintf(expectedInput, sizeof(expectedInput), "laser_fe_power");
        readAndCheckFromFile(file, expectedInput, &laser_fe_I);
        laser_fe_I *= _uW/_cm2;            // Have to get the units right!

        snprintf(expectedInput, sizeof(expectedInput), "laser_fe_s3");
        readAndCheckFromFile(file, expectedInput, &laser_fe_s3_over_s0);

        snprintf(expectedInput, sizeof(expectedInput), "laser_fe_linewidth");
        readAndCheckFromFile(file, expectedInput, &laser_fe_linewidth);
        laser_fe_linewidth *= _MHz;            // Have to get the units right!

        snprintf(expectedInput, sizeof(expectedInput), "laser_fe_nomTune");
        readAndCheckFromFile(file, expectedInput, &nominalSublevelTune2_fe);

        snprintf(expectedInput, sizeof(expectedInput), "laser_fe_detune");
        readAndCheckFromFile(file, expectedInput, &laser_fe_detune);
        laser_fe_detune *= _MHz;            // Have to get the units right!

        snprintf(expectedInput, sizeof(expectedInput), "laser_fe_offTime");
        readAndCheckFromFile(file, expectedInput, &laser_fe_offTime);
        laser_fe_offTime *= _ns;            // Have to get the units right!

        snprintf(expectedInput, sizeof(expectedInput), "laser_ge_power");
        readAndCheckFromFile(file, expectedInput, &laser_ge_I);
        laser_ge_I *= _uW/_cm2;            // Have to get the units right!

        snprintf(expectedInput, sizeof(expectedInput), "laser_ge_s3");
        readAndCheckFromFile(file, expectedInput, &laser_ge_s3_over_s0);

        snprintf(expectedInput, sizeof(expectedInput), "laser_ge_linewidth");
        readAndCheckFromFile(file, expectedInput, &laser_ge_linewidth);
        laser_ge_linewidth *= _MHz;            // Have to get the units right!

        snprintf(expectedInput, sizeof(expectedInput), "laser_ge_nomTune");
        readAndCheckFromFile(file, expectedInput, &nominalSublevelTune2_ge);

        snprintf(expectedInput, sizeof(expectedInput), "laser_ge_detune");
        readAndCheckFromFile(file, expectedInput, &laser_ge_detune);
        laser_ge_detune *= _MHz;            // Have to get the units right!

        snprintf(expectedInput, sizeof(expectedInput), "laser_ge_offTime");
        readAndCheckFromFile(file, expectedInput, &laser_ge_offTime);
        laser_ge_offTime *= _ns;            // Have to get the units right!

        snprintf(expectedInput, sizeof(expectedInput), "B_z");
        readAndCheckFromFile(file, expectedInput, &B_z);
        B_z *= _G;            // Have to get the units right!

        snprintf(expectedInput, sizeof(expectedInput), "B_x");
        readAndCheckFromFile(file, expectedInput, &B_x);
        B_x *= _G;            // Have to get the units right!

        snprintf(expectedInput, sizeof(expectedInput), "tilt");
        readAndCheckFromFile(file, expectedInput, &tilt);
        //        printf("Method = %s\n", method.c_str());
      } else {                  // File does not exist
        printf("File %s does not exist\n", inFile);
        exit(1);
      }
    } else {
      isotope = argv[1];
      if (argc > 2) Je2 = atoi(argv[2]);
      if (argc > 3) tmax = atof(argv[3])*_ns;
      if (argc > 4) dt = atof(argv[4])*_ns;
      if (argc > 5) B_z = atof(argv[5])*_G;
      if (argc > 6) B_x = atof(argv[6])*_G;
      if (argc > 7) laser_fe_detune = atof(argv[7]) *_MHz;
      if (argc > 8) laser_ge_detune = atof(argv[8]) * _MHz;
      if (argc > 9) {
        laser_fe_linewidth = atof(argv[9]) *_MHz;
        laser_ge_linewidth = atof(argv[9]) *_MHz;
      }
      if (argc > 10) laser_fe_I = atof(argv[10]) * _uW/_cm2;
      if (argc > 11) laser_ge_I = atof(argv[11]) * _uW/_cm2;
      if (argc > 12) laser_fe_s3_over_s0 = atof(argv[12]);
      if (argc > 13) laser_ge_s3_over_s0 = atof(argv[13]);
      if (argc > 14) snprintf(outFile, sizeof(outFile), "%s", argv[14]);
      if (argc > 15) method = argv[15];
      if (argc > 16) tilt = atof(argv[16]);
    }
  }
  OpticalPumping pumper;

  // printf("*********************************************\n");
  // printf("Warning...fixing laser_ge/laser_fe = 0.5\n");
  // printf("*********************************************\n");
  // laser_ge_I = laser_fe_I / 2.0;
  int status = pumper.pump(isotope, method, tmax, dt, zCoherences,
                           hfCoherences_ex, hfCoherences_gr, Je2,
                           nominalSublevelTune2_fe, nominalSublevelTune2_ge,
                           laser_fe_I, laser_ge_I, laser_fe_detune,
                           laser_ge_detune, laser_fe_linewidth,
                           laser_ge_linewidth, laser_fe_s3_over_s0,
                           laser_ge_s3_over_s0, laser_fe_offTime,
                           laser_ge_offTime, B_z, B_x, outFile, tilt);
  if (!op_batch) printf("\nCompleted with status = %d\n\n", status);
  return status;
  }
