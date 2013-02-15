// Authors: Benjamin Fenker 2013

#include <cstdio>
#include <vector>

#include "Math/MinimizerOptions.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TFrame.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TList.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVirtualFitter.h"

#define kNumLines 20000
using std::vector;
using std::pow;
using std::snprintf;

// Bool_t short = kTRUE;                    // false --> long
Double_t lastComputedPower = 0.0, lastComputedS3 = 0.0, lastComputedDetune = 0.0;
Double_t eps = 1E-6;
vector<Double_t> esPopModel;
vector<Double_t> timeModel;
vector<Double_t> polarization;
vector<Double_t> alignment;
UInt_t searchI;
Int_t runNumber, bgRunNumber;
Double_t negL;
Double_t posL;
Double_t posR;


Double_t op_tmax;



char *dataFilename = "";

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

Int_t ii = 0;
Int_t ndf = 0;
Double_t dTime[kNumLines], dFlou[kNumLines];
Double_t bgTime[kNumLines], bgFlou[kNumLines];
Double_t timeErr[kNumLines], flouErr[kNumLines];

Float_t opModel(Double_t time, Double_t *par) {
  numFuncCalls++;
  // printf("%g\n", par[0]);
  Double_t power, s3, norm, background, timeCorr, detune;
  Double_t y;
  power = par[kPower];
  s3 = par[kS3];
  norm = par[kNormalization];
  background = par[kBackground];
  timeCorr = time + par[kt_0];       // Corrected time
  // printf("Corrected t = %12.10G to %12.10G with offset of %12.10G\n", time, timeCorr, par[kt_0]);
  detune = par[kDetune];

  // Penalties work better than bounds on parameters
  if (norm < 0.0) return -20.0;         // Big penalty
  if (s3 <= 0.0 || s3 > 1.0) return -25.0; // Big penalty
  if (power < 0.0) return -30.0;          // Big penalty

  if (timeCorr > 0) {
    // printf("At time %g...", timeCorr);
    // if ((fabs(power - lastComputedPower)) > eps ||
    //     (fabs(s3 - lastComputedS3)) > eps ||
    //     (fabs(detune - lastComputedDetune))) {
    if ((power != lastComputedPower) ||
        (s3 != lastComputedS3) ||
        (detune != lastComputedDetune))  {
      numOPCodeCalls++;
      // Call my cpp-code again
      char command[500];
      sprintf(command,
               "./opticalPumping 41K 1 %g 0.1 2.0 %g %g 0.2 %g %g %g %g %s",
	      pow(10, 9)*op_tmax, detune, detune/2, power, power/2, -s3, -s3,
	      dataFilename);
      // 5 ns time steps work for RE, 0.1 for OBE
      // 2.2e+05 for tmax works for shorter test data sets
      printf("%s\n", command);
      gSystem -> Exec(command);

      FILE *fitFile = fopen(dataFilename, "r");
      Float_t ttemp, ftemp, ptemp, atemp;
      timeModel.clear();
      esPopModel.clear();
      polarization.clear();
      while (!feof(fitFile)) {
        fscanf(fitFile, " %g %g %g %g\n", &ttemp, &ftemp, &ptemp, &atemp);
        timeModel.push_back(static_cast<Double_t>(ttemp));
        esPopModel.push_back(static_cast<Double_t>(ftemp));
        polarization.push_back(static_cast<Double_t>(ptemp));
        alignment.push_back(static_cast<Double_t>(atemp));
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
      if (timeDiff > lastTimeDiff) {    // I expect the time difference to start
        // high and slowly get lower as it approaches zero.  Since I'm doing
        // everything in abs it will get near zero and start to go back up.
        // Both of these processes will be monotonic
        // printf("For index %d, I found time = %12.10G while looking for %12.10G\n", searchI,
        //        timeModel[searchI], timeCorr);
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
    y = ((norm*esPopModel[searchI])+background) + backgroundPlus;
  } else {                              // timeCorr < 0
    y = background + backgroundMinus;
    // just the background
  }
  // printf("Returning %12.10G\n", y);
  return y;
}

void calc_chi_square(Int_t &npar, Double_t *gin, Double_t &f,
                            Double_t *par, Int_t iflag) {
  double chisq = 0.0;
  ndf = 0;
  for (Int_t j = 0; j < ii; j++) {
    if (dTime[j] >= -2.0E-05 && dTime[j] < 2.2E-04) {
      // printf("For j = %d, time = %8.6G\t Comparing %8.6G\t to %8.6G with error %8.6G\n",
      //        j, dTime[j], dFlou[j], opModel(dTime[j], par), flouErr[j]);
      Double_t delta = (dFlou[j] - opModel(dTime[j], par))/flouErr[j];
      chisq += delta*delta;
      ndf++;
    }
  }
  // printf("power = %11.9G, s3 = %11.9G, norm = %11.9G, t_0 = %11.9G, detune = %11.9G, chi2 = %11.9G\n",
  //        par[0], par[1], par[2], par[kt_0], par[kDetune], chisq);
  f = chisq;
  return;
}

Float_t opFit(Double_t *time, Double_t *par) {
  return opModel(*time, par);
}

void fitOP_combined() {
  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(11);
  cout << "Enter run number.  Options are 99, 100, 105, 106" << endl;
  cin >> runNumber;
  switch (runNumber) {
  case 99:
    bgRunNumber = 101; negL = -0.1E-3; posL = 0.04E-3; posR = 0.0011;
    break;
  case 100:
    bgRunNumber = 102; negL = -0.1E-3; posL = 0.04E-3; posR = 0.0011;
    break;
  case 105:
    bgRunNumber = 103; negL = -20E-6; posL = 0.04E-3; posR = 0.22E-3;
    break;
  case 106:
    bgRunNumber = 104; negL = -20E-6; posL = 0.04E-3; posR = 0.22E-3;
    break;
  default:
    cout << "You must choose a valid run number." << endl;
    exit(1);
  }

  gStyle -> SetOptFit(1111);

  Bool_t shortie = kFALSE;                    // false --> long


  char fileName[100];
//   snprintf(fileName, sizeof(fileName), "testData/NewFile%d.dat", runNumber);
  sprintf(fileName, "testData/NewFile%d.dat", runNumber);
  FILE *dataFile = fopen(fileName, "r");
//   snprintf(fileName, sizeof(fileName), "testData/NewFile%d.dat", bgRunNumber);
  sprintf(fileName, "testData/NewFile%d.dat", bgRunNumber);
  printf("%s\n", fileName);
  FILE *bgFile = fopen(fileName, "r");
  Float_t ttemp, ftemp;
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

  TGraphErrors *dataGoodErrs = new TGraphErrors(ii, dTime, dFlou,
                                                timeErr, flouErr);
  char *tit = "";
  sprintf(tit, "Fluorescence-Run%d", runNumber);
  dataGoodErrs -> SetTitle(tit);     // Shows up on the canvas
  dataGoodErrs -> SetName(tit);      // Shows up in the fit panel
  dataGoodErrs -> Draw("A*");
  // Set up fitting function.....
  if (shortie) {
    op_tmax = (3 * posL);
  } else {
    op_tmax = posR;
  }
  printf("OP tmax = %g\n", op_tmax);

  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 6); // 6 = max params
  fitter -> SetParameter(0, "power"        , 421.0, 1.0  , 0   , 0);
  fitter -> SetParameter(1, "s3"           , 0.99 , 0.001, 0   , 1.0);
  fitter -> SetParameter(2, "normalization", 0.130, 0.001, 0   , 100.0);
  fitter -> SetParameter(3, "vert offset  ", 0.0  , 0.001, 0   , 0);
  if (runNumber == 99 || runNumber == 100) {
    fitter -> SetParameter(4, "t0"           , 0.0  , 0.001, -1.0, 1.0);
    fitter -> FixParameter(kt_0);
  } else {
    fitter -> SetParameter(4, "t0"           ,-1E-6  , 0.001, -1.0, 1.0);
  }
  fitter -> SetParameter(5, "detune"       , 2.000, 0.001, 0   , 0);

  // fitter -> FixParameter(kS3);
  // fitter -> FixParameter(kPower);
  fitter -> FixParameter(kBackground);
  fitter -> FixParameter(kDetune);
  // fitter -> FixParameter(kNormalization);

  fitter -> SetFCN(calc_chi_square);

  // Set tolerances and precision stuff
  // ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(10);
  // ROOT::Math::MinimizerOptions::SetDefaultPrecision(1E-9);
  // ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10E30);
  ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(100000);

  sprintf(dataFilename, "opData%d.dat", runNumber);

  double arglist[10];

  arglist[0] = 2;                       // print level
  fitter -> ExecuteCommand("SET PRINT", arglist, 2);
  
  arglist[0] = 100000;                    // Max function calls
  arglist[1] = 1;                       // Tolerance
  fitter -> SetMaxIterations(100000);
  printf("%d\n", fitter -> GetMaxIterations());
  fitter -> ExecuteCommand("MIGRAD", arglist, 2);

//   arglist[0] = 100000;
//   arglist[1] = kS3+1;
//   fitter -> ExecuteCommand("MINOS", arglist, 3);
  //  get result
  double minParams[6];
  double parErrors[6];
  for (int i = 0; i < 6; ++i) {  
    minParams[i] = fitter -> GetParameter(i);
    parErrors[i] = fitter -> GetParError(i);
  }
  double chi2, edm, errdef; 
  int nvpar, nparx;
  fitter -> GetStats(chi2, edm, errdef, nvpar, nparx);
  ndf -= nvpar;

  TF1 *fitFunc = new TF1("opFit", opFit, negL, op_tmax, 6);
  for (int i = 0; i < 6; i ++) {
    fitFunc -> SetParName(i, fitter -> GetParName(i));
  }

  fitFunc -> SetChisquare(chi2);
  fitFunc -> SetNDF(ndf);
  fitFunc -> SetParameters(minParams);
  fitFunc -> SetParErrors(parErrors);
  fitFunc -> FixParameter(kBackground, fitFunc -> GetParameter(kBackground));
  fitFunc -> FixParameter(kDetune, fitFunc -> GetParameter(kDetune));
  fitFunc -> FixParameter(kt_0, fitFunc -> GetParameter(kt_0));
  fitFunc -> SetLineColor(2);
  dataGoodErrs -> GetListOfFunctions() -> Add(fitFunc);

  /*
    Some crazy stuff I did to draw them on the same canvas unscaled
  TPad *pad = new TPad("pad","",0,0,1,1);
  pad->SetFillColor(42);
  pad->SetGrid();
  pad->Draw();
  pad->cd();

  TH1F *hr = c1->DrawFrame(dTime[0],-0.1,op_tmax,1.0);
  hr->SetXTitle("Time [s]");
  hr->SetYTitle("Polarization [%]");
  pad->GetFrame()->SetFillColor(21);
  pad->GetFrame()->SetBorderSize(12);


  polGraph -> Draw("LP");

  polCan -> cd();
  TPad *overlay = new TPad("overlay","",0,0,1,1);
  overlay->SetFillStyle(4000);
  overlay->SetFillColor(0);
  overlay->SetFrameFillStyle(4000);
  overlay->Draw();
  overlay->cd();

  Double_t xmin = pad->GetUxmin();
  Double_t ymin = 0.0;
  Double_t xmax = pad->GetUxmax();
  Double_t ymax = 1.0;
  TH1F *hframe = overlay->DrawFrame(xmin,ymin,xmax,ymax);
  hframe -> GetYaxis() -> SetTitle("Excited State Population [%]");
  hframe -> GetXaxis() -> SetLabelOffset(99);
  hframe -> GetYaxis() -> SetLabelOffset(99);
  
  fitFunc -> Draw();
  TGaxis *axis = new TGaxis(xmax,ymin,xmax, ymax,ymin,ymax,510,"+L");
  axis->SetLineColor(kRed);
  axis->SetLabelColor(kRed);
  axis->Draw();
  */

   for (int i = 0; i < 6; i++) {
    printf("%15s : %12.10G +/- %12.10G", fitFunc -> GetParName(i),
           fitFunc -> GetParameter(i), fitFunc -> GetParError(i));
    printf("\n");
  }
  cout << "chi2/ndf = " << chi2 << "/" << ndf << " ---> prob = ";
  cout << fitFunc -> GetProb() << endl << endl;
  
  Double_t res[kNumLines], dres[kNumLines];
  for (int i = 0; i < ii; i++) {
    res[i] = dFlou[i] - opModel(dTime[i], minParams);
    dres[ii] = dFlou[ii];
  }
  
  TGraphErrors *residualsTGE = new TGraphErrors(ii, dTime, res,
                                                timeErr, flouErr);
  c1 -> cd(2);
  residualsTGE -> SetTitle("Residuals");
  residualsTGE -> SetName("Residuals");
  residualsTGE -> Draw("A*");

  TCanvas *polCan = new TCanvas("polarization", "pol", 40, 30, 700, 700);
  polCan -> cd();
  
  Double_t flouScaled[kNumLines];
  Double_t dflouScaled[kNumLines];
  Double_t fitFuncScaled[kNumLines];
  for (int i = 0; i < ii; i++) {
    flouScaled[i] = dFlou[i]*333;
    dflouScaled[i] = flouErr[i]*333;
    fitFuncScaled[i] = (fitFunc -> Eval(dTime[i]))*333;
  }
  Double_t pol[kNumLines], polTime[kNumLines], polScaled[kNumLines];
  Double_t ali[kNumLines], aliScaled[kNumLines];
  sprintf(dataFilename, "opData%d.dat", runNumber);
  FILE *fitFile = fopen(dataFilename, "r");
  rewind(fitFile);
  Float_t ptemp, atemp;
  Int_t z = 0;
  while (!feof(fitFile)) {
    fscanf(fitFile, " %g %g %g %g\n", &ttemp, &ftemp, &ptemp, &atemp);
    pol[z] = fabs(ptemp);
    ali[z] = fabs(atemp);
    polScaled[z] = fabs(ptemp/333.0);
    aliScaled[z] = fabs(atemp/333.0);
    polTime[z] = pow(10, -6)*ttemp;

    // printf("z = %d, polTime = %8.6G, pol = %8.6G\n", z, polTime[z], pol[z]);
    z++;
  }

  TGraphErrors *dataGoodErrsScaled = new TGraphErrors(ii, dTime, flouScaled,
                                                      timeErr, dflouScaled);
  TGraph *fitFuncScaledGraph = new TGraph(ii, dTime, fitFuncScaled);
  TGraph *polGraph = new TGraph(z, polTime, pol);
  TGraph *aliGraph = new TGraph(z, polTime, ali);
  dataGoodErrsScaled -> GetXaxis() -> SetTitle("Time [s]");
  dataGoodErrsScaled -> GetYaxis() ->
    SetTitle("Fluorescence or Polarizatin [%]");
  dataGoodErrsScaled -> GetYaxis() -> SetLabelSize(0.03);
  sprintf(tit, "Fluorescence-Run%d", runNumber);
  dataGoodErrsScaled -> SetTitle(tit);
  dataGoodErrsScaled -> Draw("A*");
  fitFuncScaledGraph -> SetLineColor(2);
  fitFuncScaledGraph -> SetTitle("Fit");
  fitFuncScaledGraph -> Draw("SAME");
  polGraph -> SetLineColor(4);
  polGraph -> SetLineStyle(2);
  polGraph -> SetLineWidth(2);
  polGraph -> SetTitle("Nuclear polarization");
  polGraph -> Draw("SAME");
  aliGraph -> SetLineColor(6);
  aliGraph -> SetLineStyle(3);
  aliGraph -> SetLineWidth(2);
  aliGraph -> SetTitle("Nuclear alignment");
  aliGraph -> Draw("SAME");

  char *outfitName = "";
  sprintf(outfitName, "funcAndPol%d.pdf", runNumber);
  polCan -> SetGridy(1);
  polGraph -> GetXaxis() -> SetRangeUser(-0.005e-03, 0.1e-03);
  fitFuncScaledGraph -> GetXaxis() -> SetRangeUser(-0.005e-03, 0.1e-03);
  dataGoodErrsScaled -> GetXaxis() -> SetRangeUser(-0.005e-03, 0.1e-03);
  dataGoodErrsScaled -> GetYaxis() -> SetRangeUser(0.0, 1.01);
//   polCan -> BuildLegend(0.5, 0.5, 0.88, 0.71);
  TLegend *leg = new TLegend(0.60, 0.58, 0.96, 0.84);
  leg -> AddEntry(dataGoodErrsScaled, "Data", "PE");
  leg -> AddEntry(fitFuncScaledGraph, "Fit" , "L");
  leg -> AddEntry(polGraph, "P = #frac{<I_{z}>}{I}", "L");
  leg -> AddEntry(aliGraph,
	  "Align = #frac{I(I+1)-3<I_{z}^{2}>}{I(2I-1)}", "L");
  leg -> SetTextSize(0.03);
  leg -> Draw();
  polCan -> SaveAs(outfitName);

  c1 -> cd(1);
  TGraph *polGraphScaled = new TGraph(z, polTime, polScaled);
  TGraph *aliGraphScaled = new TGraph(z, polTime, aliScaled);
  polGraphScaled -> SetLineColor(4);
  polGraphScaled -> SetLineStyle(2);
  aliGraphScaled -> SetLineColor(6);
  aliGraphScaled -> SetLineStyle(3);
  
  polGraphScaled -> Draw("SAME");
  sprintf(outfitName, "fit%d.pdf", runNumber);
  c1 -> SaveAs(outfitName, "pdf");
}
