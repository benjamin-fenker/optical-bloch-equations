// Authors: Benjamin Fenker 2013

#include <vector>

#include "Math/MinimizerOptions.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TSystem.h"

#define kNumLines 20000
using std::vector;
using std::pow;

// Bool_t short = kTRUE;                    // false --> long
Double_t lastComputedPower = 0.0, lastComputedS3 = 0.0, lastComputedDetune = 0.0;
Double_t eps = 1E-6;
vector<Double_t> esPopModel;
vector<Double_t> timeModel;
UInt_t searchI;
Int_t runNumber, bgRunNumber;
Double_t negL;
Double_t posL;
Double_t posR;


Double_t op_tmax;

// Fit parameter numbers
const Int_t kPower = 0;
const Int_t kS3 = 1;
const Int_t kNormalization = 2;
const Int_t kBackground = 3;
const Int_t kt_0 = 4;
const Int_t kDetune = 5;

Double_t backgroundMinus = 0.0;
Double_t backgroundPlus = 0.0;

Int_t numFuncCalls = 0;
Int_t numOPCodeCalls = 0;

Float_t opModel(Double_t *time, Double_t *par) {
  numFuncCalls++;
  // printf("%g\n", par[0]);
  Double_t power, s3, norm, background, timeCorr, detune;
  power = par[kPower];
  s3 = par[kS3];
  norm = par[kNormalization];
  background = par[kBackground];
  timeCorr = time[0] + par[kt_0];       // Corrected time
  detune = par[kDetune];

  // Penalties work better than bounds on parameters
  if (norm < 0.0) return -20.0;         // Big penalty
  if (s3 < 0.0 || s3 > 1.0) return -25.0; // Big penalty
  if (power < 0.0) return -30.0;          // Big penalty

  if (timeCorr > 0) {
    // printf("At time %g...", timeCorr);
    if ((fabs(power - lastComputedPower)) > eps ||
        (fabs(s3 - lastComputedS3)) > eps ||
        (fabs(detune - lastComputedDetune))) {
      numOPCodeCalls++;
      // Call my cpp-code again
      char command[200];
      snprintf(command, sizeof(command),
               "./opticalPumping 41K 1 %g 5 2.0 %g %g 0.2 %g %g %g %g",
               op_tmax*pow(10, 9), detune, detune, power, power, s3, s3);
      printf("%s\n", command);
      gSystem -> Exec(command);

      FILE *fitFile = fopen("opData.dat", "r");
      Float_t ttemp, ftemp;
      timeModel.clear();
      esPopModel.clear();
      while (!feof(fitFile)) {
        fscanf(fitFile, " %g %g\n", &ttemp, &ftemp);
        timeModel.push_back(static_cast<Double_t>(ttemp));
        esPopModel.push_back(static_cast<Double_t>(ftemp));
      }
      fclose(fitFile);
      lastComputedPower = power;
      lastComputedS3 = s3;
      lastComputedDetune = detune;
    }

    Bool_t found = false;
    Double_t lastTimeDiff = pow(10, 6);   // Start with something huge
    searchI = 0;
    // printf("Looking for time = %6.4G...and found ", pow(10,6) *timeCorr);
    while (!found && (searchI < timeModel.size())) {
      Double_t timeM = pow(10, -6) * timeModel[searchI];
      Double_t timeDiff = fabs(timeM - timeCorr);
      // Setting the tolerance here is tricky
      // if (fabs(timeM - timeCorr) < pow(10, -6)) {
      if (timeDiff > lastTimeDiff) {    // I expect the time difference to start
        // high and slowly get lower as it approaches zero.  Since I'm doing
        // everything in abs it will get near zero and start to go back up.
        // Both of these processes will be monotonic
        found = true;
        searchI--;                      // Because we've overshot by one
      } else {
        lastTimeDiff = timeDiff;
        searchI++;
      }
    }
    // searchI--;
    // printf("Returning %g\n", esPopModel[searchI]);
    // printf("time = %6.4G\n", timeModel[searchI]);
    return ((norm*esPopModel[searchI])+background) + backgroundPlus;
  } else {                              // timeCorr < 0
    return background + backgroundMinus;
    // just the background
  }
}

void fitOP() {

  runNumber = 99;  bgRunNumber = 101; negL = -0.1E-3; posL = 0.04E-3; posR = 0.0011;
  // runNumber = 100;  bgRunNumber = 102; negL = -0.1E-3; posL = 0.04E-3; posR = 0.0011;
  // runNumber = 105;  bgRunNumber = 103; negL = -20E-6; posL = 0.04E-3; posR = 0.22E-3;
  // runNumber = 106;  bgRunNumber = 104; negL = -20E-6; posL = 0.04E-3; posR = 0.22E-3;


  Bool_t shortie = kTRUE;                    // false --> long


  char fileName[100];
  snprintf(fileName, sizeof(fileName), "testData/NewFile%d.dat", runNumber);
  FILE *dataFile = fopen(fileName, "r");
  snprintf(fileName, sizeof(fileName), "testData/NewFile%d.dat", bgRunNumber);
  printf("%s\n", fileName);
  FILE *bgFile = fopen(fileName, "r");
  Double_t dTime[kNumLines], dFlou[kNumLines];
  Double_t bgTime[kNumLines], bgFlou[kNumLines];
  Double_t timeErr[kNumLines], flouErr[kNumLines];
  Float_t ttemp, ftemp;
  Int_t ii = 0;
  TH1D *residualsHist = new TH1D("res", "residues", 5, -0.00025, 0.00025);
  while (!feof(dataFile)) {
    fscanf(dataFile, " %g %g\n", &ttemp, &ftemp);
    // printf("Read %g and %g\n", ttemp, ftemp);
    dTime[ii] = static_cast<Double_t>(ttemp);
    dFlou[ii] = static_cast<Double_t>(ftemp);
    timeErr[ii] = 0.0;                  // No appreciable error on time
    flouErr[ii] = 7.4E-5;
    residualsHist -> Fill(dFlou[ii] - 0.0008499);
    // Result of fit to flat background

    fscanf(bgFile, " %g %g\n", &ttemp, &ftemp);
    // printf("Read %g and %g\n", ttemp, ftemp);
    bgTime[ii] = static_cast<Double_t>(ttemp);
    bgFlou[ii] = static_cast<Double_t>(ftemp);
    ii++;
  }

  // Fit the background data in two regions.  t < 0 and t > 0
  Bool_t plotBackground = false;
  TGraph *bg = new TGraph(ii, bgTime, bgFlou);
  if (plotBackground) {
    TCanvas *bgCan = new TCanvas("background", "bg", 40, 30, 700, 700);
    bgCan -> cd();
    bg -> Draw("A*");
  }
  TF1 *negFit = new TF1("negFit", "[0]", negL, 0.0);
  bg -> Fit(negFit, "R");
  backgroundMinus = negFit -> GetParameter(0);
  TF1 *posFit = new TF1("posFit", "[0]", posL, posR);
  bg -> Fit(posFit, "R");
  backgroundPlus = posFit -> GetParameter(0);

  TCanvas *c1 = new TCanvas("fit results", "fit", 40, 30, 700, 700);
  c1 -> Divide(1, 2);
  c1 -> cd(1);

  // residualsHist -> Draw("HIST");
  TGraphErrors *data = new TGraphErrors(ii, dTime, dFlou, timeErr, flouErr);
  // Fit the scatter of the data points in the long flat tail region
  TF1 *lin = new TF1("scatter", "[0]", posL, posR);
  TFitResultPtr p = data -> Fit(lin, "SRQ");
  Double_t eAdjust = sqrt((p -> Chi2()) / (p -> Ndf()));
  printf("Adjusting errors by %g\n", eAdjust);
  for (Int_t jj = 0; jj < ii; jj++) {
    flouErr[jj] *= eAdjust;
  }

  
  // Set tolerances and precision stuff
  ROOT::Math::MinimizerOptions::SetDefaultPrecision(1E-9);
  // ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10E30);
  // ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(10E30);


  TGraphErrors *dataGoodErrs = new TGraphErrors(ii, dTime, dFlou,
                                                timeErr, flouErr);
  dataGoodErrs -> SetTitle("Fluorescence");     // Shows up on the canvas
  dataGoodErrs -> SetName("Fluorescence");      // Shows up in the fit panel
  dataGoodErrs -> Draw("A*");
  // Set up fitting function.....
  if (shortie) {
    op_tmax = (3 * posL);
  } else {
    op_tmax = posR;
  }
  printf("OP tmax = %g\n", op_tmax);
  TF1 *fitFunc = new TF1("opFit", opModel, negL, op_tmax, 6); // long

  fitFunc -> SetParNames("power", "s3", "normalization",
                         "background", "t0", "detuning");
  // Set up power - don't use limits unless you have to!
  // fitFunc -> SetParLimits(kPower, 0.0, 10000.0);
  fitFunc -> FixParameter(kPower, 420.392);

  // Set up S3
  // fitFunc -> SetParLimits(kS3, 0.0, 1.0);
  fitFunc -> SetParameter(kS3, 0.999047);

  // Set up normalization
  // fitFunc -> SetParLimits(kNormalization, 0, 1000.0);
  fitFunc -> SetParameter(kNormalization, 0.128206);

  // Set up background
  fitFunc -> FixParameter(kBackground, 0.0);

  // Set up time offset
  // fitFunc -> SetParameter(kt_0, -9.0618E-7);
  // fitFunc -> FixParameter(kt_0, -9.07464E-7);
  fitFunc -> FixParameter(kt_0, 0.0);

  // Set up the detuning
  fitFunc -> SetParameter(kDetune, -0.738384);

  // ...and do the fit
  TFitResultPtr r = dataGoodErrs -> Fit(fitFunc, "SBR");
  // Int_t fitStatus = r;
  printf("Called the fit function %d times\n", numFuncCalls);
  printf("Called the OP code %d times\n", numOPCodeCalls);
  r -> Print("V");
  // Plot the residuals
  Double_t residual[kNumLines];
  for (Int_t i = 1; i < ii; i++) {
    residual[i] = dFlou[i] - (fitFunc -> Eval(dTime[i]));
  }

  // TCanvas *c2 = new TCanvas("residuals", "res", 40, 30, 700, 700);
  // c2 -> cd();
  c1 -> cd(2);
  TGraphErrors *residualsTGE = new TGraphErrors(ii, dTime, residual,
                                                timeErr, flouErr);
  residualsTGE -> SetTitle("Residuals");
  residualsTGE -> SetName("Residuals");
  residualsTGE -> Draw("A*");
}
