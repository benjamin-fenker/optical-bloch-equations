// Authors: Benjamin Fenker 2013

#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TSystem.h"

#define kNumLines 20000

Double_t lastComputedPower = 0.0, lastComputedS3 = 0.0;
Double_t eps = pow(10, -6);
vector<Double_t> esPopModel;
vector<Double_t> timeModel;
UInt_t searchI;

// Fit parameter numbers
const Int_t kPower = 0;
const Int_t kS3 = 1;
const Int_t kNormalization = 2;
const Int_t kBackground = 3;

Double_t backgroundMinus = 0.0;
Double_t backgroundPlus = 0.0;

Int_t numOPCodeCalls = 0;

Float_t opModel(Double_t *time, Double_t *par) {
  // printf("%g\n", par[0]);
  Double_t power, s3, norm, background;
  power = 4*par[kPower];
  s3 = par[kS3];
  norm = par[kNormalization];
  background = par[kBackground];

  if (time[0] > 0) {
    // printf("At time %g...", time[0]);
    if ((fabs(power - lastComputedPower)) > eps ||
        (fabs(s3 - lastComputedS3)) > eps) {
      // Call my cpp-code again
      numOPCodeCalls++;
      char command[200];
      snprintf(command, sizeof(command),
               // long
               // "./opticalPumping 41K 1 1100000 5 2.0 -1 -1 0.2 %g %g %g %g",
               // short
              "./opticalPumping 41K 1 200000 5 2.0 -1 -1 0.2 %g %g %g %g",
              power, power, s3, s3);
      // printf("%s\n", command);
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
    }

    Bool_t found = false;
    Double_t lastTimeDiff = pow(10, 6);   // Start with something huge
    searchI = 0;
    // printf("Looking for time = %6.4G...and found ", pow(10,6) *time[0]);
    while (!found && (searchI < timeModel.size())) {
      Double_t timeM = pow(10, -6) * timeModel[searchI];
      Double_t timeDiff = fabs(timeM - time[0]);
      // Setting the tolerance here is tricky
      // if (fabs(timeM - time[0]) < pow(10, -6)) {
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
  } else {                              // time[0] < 0
    return background + backgroundMinus;
    // just the background
  }
}

void fitOP() {
  // gSystem -> Exec("./opticalPumping -h");
  FILE *dataFile =fopen("testData.dat", "r");
  FILE *bgFile =fopen("testData/NewFile101_tabs.dat", "r");
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
  TCanvas *bgCan = new TCanvas("background", "bg", 40, 30, 700, 700);
  TGraph *bg = new TGraph(ii, bgTime, bgFlou);
  bgCan -> cd();
  bg -> Draw("A*");
  TF1 *negFit = new TF1("negFit", "[0]", -0.0001, 0.0);
  bg -> Fit(negFit, "R");
  backgroundMinus = negFit -> GetParameter(0);
  TF1 *posFit = new TF1("posFit", "[0]", 0.0, 0.0011);
  bg -> Fit(posFit, "R");
  backgroundPlus = posFit -> GetParameter(0);

  TCanvas *c1 = new TCanvas("fit results", "fit", 40, 30, 700, 700);
  c1 -> Divide(1, 2);
  c1 -> cd(1);

  // residualsHist -> Draw("HIST");
  TGraphErrors *data = new TGraphErrors(ii, dTime, dFlou, timeErr, flouErr);
  data -> SetTitle("Fluorescence");     // Shows up on the canvas
  data -> SetName("Fluorescence");      // Shows up in the fit panel
  data -> Draw("A*");

  // Set up fitting function.....
  // TF1 *fitFunc = new TF1("opFit", opModel, -0.0001, 0.0011, 4); // long
  TF1 *fitFunc = new TF1("opFit", opModel, -0.05E-3, 0.2E-3, 4);   // short
  fitFunc -> SetParNames("power", "s3", "normalization", "background");
  // Set up power
  fitFunc -> SetParLimits(kPower, 200.0, 10000.0);
  fitFunc -> SetParameter(kPower, 1663.0/4.0);
  // Set up S3
  fitFunc -> SetParLimits(kS3, 0.0, 1.0);
  fitFunc -> SetParameter(kS3, 0.97);
  // Set up normalization
  fitFunc -> SetParLimits(kNormalization, 0, 1000.0);
  fitFunc -> SetParameter(kNormalization, 0.1828);
  // Set up background
  fitFunc -> SetParLimits(kBackground, 0.0, 1.0);
  // fitFunc -> SetParameter(kBackground, 0.000403);
  fitFunc -> FixParameter(kBackground, 0.0);

  // ...and do the fit
  TFitResultPtr r = data -> Fit(fitFunc, "BR");
  // Int_t fitStatus = r;
  printf("Called the OP code %d time\n", numOPCodeCalls);

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
