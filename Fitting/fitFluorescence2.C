// Authors: Benjamin Fenker 2013
#include <iostream>
#include <fstream>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TClass.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TVirtualFitter.h"

#include "OPFitObj.h"

using std::cin;
using std::cout;
using std::ofstream;
using std::endl;
using std::vector;

// Until I can think of a general way to calculate these, I'm putting'em up here
const Int_t numPoints = 11;

Double_t t_left[numPoints] = {0.0, 1.0, 2.0, 4.0, 9.0, 19.0, 39.0,
                                     89.0, 199.0, 499.0, 999.0};
Double_t t_right[numPoints] = {1.0, 2.0, 4.0, 9.0, 19.0, 39.0, 89.0,
                               199.0, 499.0, 999.0, 1740.0};
TVirtualFitter *fitter = 0;
static Int_t kPower = 0;
static Int_t ks3 = 1;
static Int_t kNorm = 2;
static Int_t kDetune = 3;
static Int_t kBkgLevel = 4;
static Int_t kBFieldZ = 5;
static Int_t kBFieldX = 6;

static Int_t kPower_sim = 0;
static Int_t ks3_plu = 1;
static Int_t ks3_min = 2;
static Int_t kNorm_plu = 3;
static Int_t kNorm_min = 4;
static Int_t kDetune_sim = 5;
static Int_t kBkgLevel_plu = 6;
static Int_t kBkgLevel_min = 7;
static Int_t kBFieldZ_sim = 8;
static Int_t kBFieldX_sim = 9;

ofstream chi2Progress("chi2Prog.txt");

// Fitting functions ***********************************************************
void printChi2Prog(Double_t  *par, Int_t npar, Double_t chi2) {
  chi2Progress.open("chi2Prog.txt", ios::app);
  for (Int_t i = 0; i < npar; i++) {
    chi2Progress << par[i] << "\t";
  }
  chi2Progress << chi2 << endl;
  chi2Progress.close();
}

Double_t integrate_trapezoid(Double_t t_start, Double_t t_stop,
                             vector<Double_t> t, vector<Double_t> fl) {
  Double_t fluor = 0.0;
  Double_t N = 0.0;
  for (UInt_t i = 0; i < fl.size()-1; i++) {
    if (t[i] >= t_start && (t_stop - t[i]) >= 0.0) {
      fluor += fl[i] + fl[i+1];
      N += 1.0;
    }
  }
  fluor = fluor / (2.0*N);
  return fluor;
}

void calc_chi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,
               Int_t iflag) {
  if (iflag == 1) cout << "Initializing..." << endl;
  if (iflag == 3) cout << "Cleaning up..." << endl;
  if (0) cout << gin[0];                // Gets rid of unused parameter warning
  Double_t chi2 = 0.0;
  TObjArray *objA =  (TObjArray *)fitter -> GetObjectFit();
  Double_t power, detune, Bz, Bx;
  Double_t s3[2], norm[2], bkg[2];
  switch (objA -> GetSize()) {
    case 1:                             // Sigma + or Sigma-
      power = par[kPower];
      s3[0] = par[ks3];
      s3[1]=  par[ks3];
      norm[0] = par[kNorm];
      norm[1] = par[kNorm];
      detune = par[kDetune];
      bkg[0] = par[kBkgLevel];
      bkg[1] = par[kBkgLevel];
      Bz = par[kBFieldZ];
      Bx = par[kBFieldX];
      break;
    case 2:                             // Simultaneous
      power = par[kPower_sim];
      s3[0] = par[ks3_plu];
      s3[1] = par[ks3_min];
      norm[0] = par[kNorm_plu];
      // norm[1] = par[kNorm_min];
      norm[1] = norm[0];                // Each set has same normalization
      detune = par[kDetune_sim];
      bkg[0] = par[kBkgLevel_plu];
      bkg[1] = par[kBkgLevel_min];
      Bz = par[kBFieldZ_sim];
      Bx = par[kBFieldX_sim];
      break;
    default:
      cout << "More than two Objs to fit to! ";
      cout << "(" << objA -> GetSize() << " to be exact)" << endl;
      exit(1);
  }
  Double_t tmp = 0.0;
  // OPFitObject *fitObj = (OPFitObject *)fitter -> GetObjectFit();

  for (Int_t i = 0; i < objA -> GetSize(); i++) {
    OPFitObject *fitObj = (OPFitObject *)objA -> At(i);

    fitObj -> setCurrentPower(power);
    fitObj -> setCurrentStokes(s3[i]);
    fitObj -> setCurrentDetune(detune);
    fitObj -> setCurrentBz(Bz);
    fitObj -> setCurrentBx(Bx);
    if (fitObj -> needsARerun()) fitObj -> rerunOPModel();
    for (UInt_t j = 0; j < fitObj -> expDat.size(); j++) {
      Double_t integrated = integrate_trapezoid(t_left[j], t_right[j],
                                                fitObj -> mod_t,
                                                fitObj -> mod_e);
      integrated *= norm[i];
      integrated += bkg[i];
      tmp = (integrated - fitObj->expDat[j])/(fitObj -> dexpDat[j]);
      tmp = tmp*tmp;
      chi2 += tmp;

      tmp = (bkg[i] - fitObj->expBkg[j])/(fitObj->dexpBkg[j]);
      tmp = tmp*tmp;
      chi2 += tmp;
    }
  }
  printChi2Prog(par, npar, chi2);
  f = chi2;
}
// *****************************************************************************

// Helper functions for formatting *********************************************
void formatDataTGraph(TGraph *tg) {
  tg -> SetTitle("Fluorescence fitting");
  tg -> GetXaxis() -> SetTitle("Time [us]");
  tg -> GetYaxis() -> SetTitle("Fluorescence");
  tg -> SetMarkerStyle(4);
  tg -> SetMarkerSize(0.4);
}
// *****************************************************************************

void fitFluorescence2(Double_t sGuess, Double_t fGuess, Double_t xGuess) {
  Double_t pGuess = 700.0; // mW/cm^2
  //    Double_t sGuess = 0.965;        // Stokes (s3)
  Double_t nGuess = 40.00;        // Normalization
  Double_t dGuess = 2.000;        // Detuning
  Double_t bGuess = 0.474312;        // Background level (From linear fit)
//   Double_t fGuess = 1.900;        // B-Field_z
  //    Double_t xGuess = 0.000;        // B-Field_x (transverse)

  // Load my classes ***********************************************************
  if (!TClass::GetDict("OPFitObject")) {
    gROOT->ProcessLine(".L OPFitObj.cpp++");
  }
  // ***************************************************************************

  // Setup graphics and I/O*****************************************************
  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(11);
  chi2Progress.open("chi2Prog.txt");
  // ***************************************************************************

  // Variable declarations *****************************************************
  Int_t polChoice, lasChoice, methodChoice;           // Input choices
  vector<Double_t> t_center(numPoints, 0.0);
  vector<Double_t> dt_center(numPoints, 0.0); // Times and bin widths from JB
  // All the columns of data from JB
  vector<Double_t> sp(numPoints, 0.0), dsp(numPoints, 0.0),
      spbkg(numPoints, 0.0), dspbkg(numPoints, 0.0), sm(numPoints, 0.0),
      dsm(numPoints, 0.0), smbkg(numPoints, 0.0), dsmbkg(numPoints, 0.0);
  ifstream dataFile;
  vector<Double_t> arglist(100, 0.0);             // For MINUIT
  OPFitObject *fitObj_p, *fitObj_m;
  TCanvas *fitCan;
  TGraphErrors *dataG_p, *bkgG_p, *dataG_m, *bkgG_m;
  TObjArray *graphsToDraw, *fitObjsToSend;
  TGraph *fitG;
  vector<Double_t> fitVec;
  TF1 *fitBG;
  TPaveText *paves;
  Double_t chi2, edm, errdef;
  Int_t ndf, nvpar, nparx;
  char title[100];
  // ***************************************************************************

  // Get user input ************************************************************
  cout << "Choose polarization: enter >0 for sigma+ and <0 for sigma-" << endl;
  cout << "Enter 0 to simultaneously fit both" << endl;
  cin >> polChoice;
  if (polChoice != 0) polChoice = polChoice / abs(polChoice); // Ensures |p| = 1
  if (polChoice > 0)  cout << "You chose sigma+" << endl;
  if (polChoice < 0)  {
    cout << "You chose sigma-" << endl;
    dGuess *= -1.0;
  }
  if (polChoice == 0) cout << "You chose both" << endl;

  cout << "Choose Magnetic field state 1 or 2 (see laserScheme.pdf for what";
  cout <<" this means)." << endl;
  cin >> lasChoice;
  if (lasChoice != 1 && lasChoice != 2) {
    cout << "NO! You cannot pick " << lasChoice << endl;
    exit(1);
  }

  cout << "Choose calculation method.  RE is faster but DM does"
       << "transverse fields. 0 for RE and 1 for DM " << endl;
  cin >> methodChoice;
  if (methodChoice != 0 && methodChoice != 1) {
    cout << "Cannot choose " << methodChoice << endl;
    exit(1);
  }
  cout << "Choose which parameters to fix (others will be varied): " << endl;
  cout << "Enter -1 when done" << endl;
  if (polChoice != 0) {
    cout << "Power: " << kPower << endl << "s3: " << ks3 << endl
	 << "Normalization: " << kNorm << endl << "Detune: " << kDetune
	 << endl << "Background: " << kBkgLevel << endl << "B-field_z: "
	 << kBFieldZ << endl << "B-field_x: " << kBFieldX << endl;
  } else {
    cout << "Power (both): " << kPower_sim << endl << "s3 (+): " << ks3_plu
	 << endl << "s3 (-): " << ks3_min << endl << "Normalization (+): " << kNorm_plu
	 << endl << "Normalization (-): " << kNorm_min << endl << "Detune: " << kDetune_sim
	 << endl << "Background (+): " << kBkgLevel_plu << endl << "Background (-): "
	 << kBkgLevel_min << endl << "B-fieldz: " << kBFieldZ_sim << endl << "B-field_x: "
	 << kBFieldX_sim << endl;
 }

  Int_t fixChoice = -1;
  vector<Int_t> fix;
  cin >> fixChoice;
  while (fixChoice > 0) {
    fix.push_back(fixChoice);
    cin >> fixChoice;
  }

  cout << "Choose MIGRAD strategy (0-2)" << endl;
  Int_t migStrat;
  cin >> migStrat;

  // ***************************************************************************

  // Read data file ************************************************************
  Int_t fileChoice = -1;
  cout << "Enter 1 for old dat with systematic errors" << endl;
  cout << "Enter 2 for new (Dec-like) data with systematic errors" << endl;
  cout << "Enter 3 for old data with statistical errors only" << endl;
  cout << "Enter 4 for new (Dec-like) data with stat errors only" << endl;
  cout << "Enter 5 for new (Dec-like) data with scaled systematic errors" << endl;
  cout << "Enter 6 for old data with scaled systematic errors" << endl;
  cin >> fileChoice;
  switch (fileChoice) {
  case 1:                             // Sigma + or Sigma-
    dataFile.open("FluorescenceScalerData_plain.dat");
    break;
  case 2:
    dataFile.open("FluorescenceScalerData1286_plain.dat");
    break;
  case 3:
    dataFile.open("FluorescenceScalerData_1148-stats.dat");
    break;
  case 4:
    dataFile.open("FluorescenceScalerData_1286-stats.dat");
    break;
  case 5:
    dataFile.open("FluorescenceScalerData1286_normalized.dat");
    break;
  case 6:
    dataFile.open("FluorescenceScalerData1148_normalized.dat");
    break;
  default:
    cout << "I said enter one or two.  Follow instructions" << endl;
    exit(1);
  }
  
  if (dataFile.is_open()) {
    // cout << "Opened file." << endl;
    for (Int_t i = 0; i < numPoints; i++) {
      if (!dataFile.eof()) {
        dataFile >> t_center[i] >> sp[i] >> dsp[i] >> spbkg[i] >> dspbkg[i]
                 >> sm[i] >> dsm[i] >> smbkg[i] >> dsmbkg[i];
        dt_center[i] = (t_right[i] - t_left[i])/2.0;

        // sp_net[i] = sp[i] - spbkg[i];
        // sm_net[i] = sm[i] - smbkg[i];
        // dsp_net[i] = sqrt((dsp[i]*dsp[i])+(dspbkg[i]*dspbkg[i]));
        // dsm_net[i] = sqrt((dsm[i]*dsm[i])+(dsmbkg[i]*dsmbkg[i]));
      } else {
        cout << "Epic fail, file ended too soon" << endl;
      } // Done failing
    }   // Done for loop
  } else {
    cout << "Epic fail, couldn't open file" << endl;
  } // Done reading file
  // cout << "done reading file!" << endl;
  // ***************************************************************************
  fitCan = new TCanvas();
  fitCan -> cd();
  fitCan -> SetLogx(1);

  TVirtualFitter::SetDefaultFitter("Minuit");
  // fitObj = new OPFitObject();
  arglist[0] = 0;

  dataG_p = new TGraphErrors(t_center.size(), &t_center[0], &sp[0], 0, &dsp[0]);
  dataG_p -> SetTitle("Sigma+ data");
  dataG_p -> SetName("Sigma+ data");

  bkgG_p = new TGraphErrors(t_center.size(), &t_center[0], &spbkg[0], 0, &dspbkg[0]);
  bkgG_p -> SetTitle("Sigma+ bkg");
  bkgG_p -> SetName("Sigma+ bkg");

  dataG_m = new TGraphErrors(t_center.size(), &t_center[0], &sm[0], 0, &dsm[0]);
  dataG_m -> SetTitle("Sigma- data");
  dataG_m -> SetName("Sigma- data");

  bkgG_m = new TGraphErrors(t_center.size(), &t_center[0], &smbkg[0], 0, &dsmbkg[0]);
  bkgG_m -> SetTitle("Sigma- bkg");
  bkgG_m -> SetName("Sigma- bkg");

  fitObj_p = new OPFitObject(sp, dsp, spbkg, dspbkg, pGuess, sGuess, dGuess,
                             2.0, fGuess, xGuess, lasChoice, polChoice);
  fitObj_m = new OPFitObject(sm, dsm, smbkg, dsmbkg, pGuess, -sGuess, dGuess,
                             2.0, fGuess, xGuess, lasChoice, polChoice);
  if (methodChoice == 0) {
    fitObj_p -> setMethod_RateEquation();
    fitObj_m -> setMethod_RateEquation();
  } else {
    fitObj_p -> setMethod_DensityMatrix();
    fitObj_m -> setMethod_DensityMatrix();
  }
  graphsToDraw = new TObjArray(2);
  fitObjsToSend = new TObjArray(1);
  if (polChoice != 0) {
    fitter = TVirtualFitter::Fitter(0, 7);
    fitter -> ExecuteCommand("CLE", 0, 0);
    fitter -> SetParameter(kPower   , "power" , pGuess, 1.0, 0.0, 10000.0);
    fitter -> SetParameter(kNorm    , "norm"  , nGuess, 1.0, 0.0, 0.00000);
    fitter -> SetParameter(kDetune  , "detune", dGuess, 0.1, 0.0, 0.00000);
    fitter -> SetParameter(kBkgLevel, "BG Lev", bGuess, 0.1, 0.0, 0.00000);
    fitter -> SetParameter(kBFieldZ , "B_z"   , fGuess, 0.1, 0.0, 0.00000);
    fitter -> SetParameter(kBFieldX , "B_x"   , xGuess, 0.1, 0.0, 0.17330);  // G (5 degree tilt)


    if (polChoice > 0) {
      fitter -> SetParameter(ks3, "s3"  , sGuess, 0.001, 0.0, 1.00000);
      graphsToDraw -> Add(dataG_p);
      graphsToDraw -> Add(bkgG_p);
      fitObjsToSend -> Add(fitObj_p);
      sprintf(title, "Sigma+, State %d", lasChoice);
    }
    if (polChoice < 0) {
      fitter -> SetParameter(ks3, "s3"  , -sGuess, 0.001, -1.00000, 0.0);
      graphsToDraw -> Add(dataG_m);
      graphsToDraw -> Add(bkgG_m);
      fitObjsToSend -> Add(fitObj_m);
      sprintf(title, "Sigma-, State %d", lasChoice);
    }
//     fitter -> FixParameter(kDetune);
//     fitter -> FixParameter(kBFieldZ);
//     fitter -> FixParameter(kBFieldX);
//     fitter -> FixParameter(ks3);

  }
  if (polChoice == 0) {
    fitObjsToSend -> Expand(2);
    fitter = TVirtualFitter::Fitter(0, 9);
    fitter -> ExecuteCommand("CLE", 0, 0);
    fitter -> SetParameter(kPower_sim   , "power"     , pGuess  , 1.0, 0.0, 10000);
    fitter -> SetParameter(ks3_plu      , "s3 (+)"    , sGuess , 0.001, 0.0, 1.00);
    fitter -> SetParameter(ks3_min      , "s3 (-)"    , -sGuess, 0.001, -1.00, 0.0);
    fitter -> SetParameter(kNorm_plu    , "norm (+)"  , nGuess, 1.0, 0.0, 0.00);
    fitter -> SetParameter(kNorm_min    , "norm (-)"  , nGuess, 1.0, 0.0, 0.00);
    fitter -> SetParameter(kDetune_sim  , "detune"    , dGuess, 0.1, 0.0, 0.00);
    fitter -> SetParameter(kBkgLevel_plu, "BG Lev (+)", bGuess, 0.1, 0.0, 0.00);
    fitter -> SetParameter(kBkgLevel_min, "BG Lev (-)", bGuess, 0.1, 0.0, 0.00);
    fitter -> SetParameter(kBFieldZ_sim , "B_z"       , fGuess, 0.1, 0.0, 0.00);
    fitter -> SetParameter(kBFieldX_sim , "B_x"       , xGuess, 0.1, 0.0, 0.00);

//     fitter -> FixParameter(kDetune_sim);
//     fitter -> FixParameter(kBFieldZ_sim);
//     fitter -> FixParameter(kNorm_min);  // Fixed to norm+
//     fitter -> FixParameter(kBFieldX_sim);

    graphsToDraw -> Add(dataG_p);
    graphsToDraw -> Add(bkgG_p);
    graphsToDraw -> Add(dataG_m);
    graphsToDraw -> Add(bkgG_m);
    fitObjsToSend -> Add(fitObj_p);
    fitObjsToSend -> Add(fitObj_m);
    fitCan -> Divide(1, 2);
    fitCan -> cd(1);
    gPad -> SetLogx(1);
  }

  for (UInt_t i = 0; i < fix.size(); i++) {
    fitter -> FixParameter(fix[i]);
  }

  Bool_t first = kTRUE;
  for (Int_t i = 0; i < graphsToDraw -> GetSize(); i++) {
    if (graphsToDraw->At(i)->InheritsFrom("TGraph")) {
      formatDataTGraph((TGraph*)graphsToDraw->At(i));
      if (polChoice == 0 && i == 0) {
        sprintf(title, "Sigma+, State %d", lasChoice);
      }
      if (polChoice == 0 && i == 1) {
        sprintf(title, "Sigma-, State %d", lasChoice);
      }
      TGraph *t = (TGraph *)graphsToDraw -> At(i);
      t -> SetTitle(title);
      if (first) {
        graphsToDraw->At(i)->Draw("AP");
        first = kFALSE;
      } else {
        graphsToDraw->At(i)->Draw("P");
      }
    }
    if (polChoice == 0 && i == 1) {
      fitCan -> cd(2);
      gPad -> SetLogx(1);
      first = kTRUE;
    }
  }

  // Fit the data ************************************************************
  fitter -> SetObjectFit(fitObjsToSend);
  fitter -> SetFCN(calc_chi2);
  arglist[0] = migStrat;
  // fitter -> ExecuteCommand("SET PRI", &arglist[0], 1);
  fitter -> ExecuteCommand("SET STR", &arglist[0], 1);
  arglist[0] = 5000;
  arglist[1] = 1;
  int j;
  cout << "Enter 0 to skip fitting, 1 to fit but skip HESSE or anything else"
       << " for both" << endl;
  cin >> j;
  if (j != 0) {
    fitter -> ExecuteCommand("MIGRAD", &arglist[0], 2);
    if (j != 1) {
      fitter -> ExecuteCommand("HESSE", &arglist[0], 1);
    }
  }
  arglist[1] = 0;
//   fitter -> ExecuteCommand("SHO COV", &arglist[0], 0);
//   fitter -> ExecuteCommand("MINOS", &arglist[0], 1);
  ndf = sm.size() + smbkg.size() - fitter -> GetNumberFreeParameters();
  if (polChoice == 0) ndf += sm.size() + smbkg.size();
  fitter -> GetStats(chi2, edm, errdef, nvpar, nparx);

  cout << "Chi2/NDF = " << chi2 << " / " << ndf << " = " << chi2/ndf
       << " --> prob = " <<   TMath::Prob(chi2, ndf) << endl;
  for (Int_t i = 0; i < fitter -> GetNumberTotalParameters(); i++) {
    cout << fitter -> GetParName(i) << ": " << fitter -> GetParameter(i)
         << " +/- " << fitter -> GetParError(i) << endl;
  }

  cout << "****************************************************" << endl;
  Double_t scaleFactor = sqrt(chi2/ndf);
  cout << "Assuming statistical uncertainties not in the model and that errors "
       << "have been overestimated." << endl << "Scaling errors by: "
       << sqrt(chi2/ndf) << " leads to: " << endl;
  for (Int_t i = 0; i < fitter -> GetNumberTotalParameters(); i++) {
    cout << fitter -> GetParName(i) << ": "
         << fitter -> GetParameter(i) << " +/- "
         << fitter -> GetParError(i)*scaleFactor << endl;
  }

  printf("Chi2/NDF = %8.6G/%d = %G --> prob = %G\n", chi2, ndf,
	 chi2/(Double_t)ndf, TMath::Prob(chi2, ndf));
  for (Int_t i = 0; i < fitter -> GetNumberTotalParameters(); i++) {
    printf("%8s: %8.6G +/- %8.6G (%8.6G)\n", fitter -> GetParName(i),
	   fitter -> GetParameter(i), fitter -> GetParError(i),
	   fitter -> GetParError(i)*scaleFactor);
  }
  // Draw the result of the fit **********************************************
  if(polChoice == 0) fitCan -> cd(1);
  vector<Double_t> norm(2, 0.0);
  vector<Double_t> bkg(2, 0.0);
  cout << "Size is " << fitObjsToSend -> GetSize() << endl;
  if (fitObjsToSend -> GetSize() == 1) {
    norm[0] = fitter -> GetParameter(kNorm);
    norm[1] = norm[0];
    bkg[0] = fitter -> GetParameter(kBkgLevel);
    bkg[1] = bkg[0];
  } else if (fitObjsToSend -> GetSize() == 2) {
    norm[0] = fitter -> GetParameter(kNorm_plu);
    norm[1] = fitter -> GetParameter(kNorm_min);
    bkg[0] = fitter -> GetParameter(kBkgLevel_plu);
    bkg[1] = fitter -> GetParameter(kBkgLevel_min);
  } else {
    cout << "Length is " << fitObjsToSend -> GetSize()
         << "  be 1 or 2" << endl;
    exit(1);
  }
  for (Int_t i = 0; i < fitObjsToSend -> GetSize(); i++) {
    Double_t b = bkg[i];
    Double_t n = norm[i];
    OPFitObject *o = (OPFitObject *)fitObjsToSend -> At(i);
    fitVec = o -> applyNormAndBackground(n, b);
    fitG = new TGraph(o -> mod_t.size(), &(o -> mod_t[0]),
                      &fitVec[0]);
    fitG -> SetLineColor(2);
    fitG -> Draw("L");
    fitBG = new TF1("bgResult", "[0]", 0.0, t_right[numPoints-1]);
    fitBG -> SetParameter(0, b);
    fitBG -> SetLineWidth(1);
    fitBG -> Draw("SAME");
    if (polChoice == 0) fitCan -> cd(2);
  }

  // Print the best fit parameters on the Canvas *****************************
  char label[100];
  paves = new TPaveText(0.5, 0.5, 0.9, 0.9, "brNDC");
  // paves = new TPaveText(1.21997, 0.938843, 3.14517, 1.45549);
  sprintf(label, "chi2/ndf =  %g/%d --> %g\n", chi2, ndf,
          TMath::Prob(chi2, ndf));
  paves -> AddText(label);
  cout << "Tots: " << fitter -> GetNumberTotalParameters() << endl;
  for (Int_t i = 0; i < fitter -> GetNumberTotalParameters(); i++) {
//     if (!fitter -> IsFixed(i)) {
      sprintf(label, "%8s: %8.6g +/- %8.6g\n", fitter -> GetParName(i),
              fitter -> GetParameter(i), fitter -> GetParError(i));
      paves -> AddText(label);
//     }
  }
  paves -> Draw();
  
  sprintf(title, "out%g.C", chi2);
  fitCan -> SaveAs(title);
  char pChar;
  if (polChoice == 1) pChar = 'p';
  if (polChoice == -1) pChar = 'm';
  sprintf(title, "outB%g-s%c-l%d.ps", (fitter -> GetParameter(kBFieldX))*100.0,
	  pChar, lasChoice);
  fitCan -> SaveAs(title);
}   // End main
