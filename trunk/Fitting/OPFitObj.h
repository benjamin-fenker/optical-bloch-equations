/* Authors: Benjamin Fenker 2013 */
#include <time.h>

#include <iostream>
#include <fstream>
#include <vector>

#include "TObject.h"
#include "TSystem.h"

using std::vector;
using std::ofstream;

class OPFitObject : public TObject {
public:
  OPFitObject();
  OPFitObject(vector<Double_t> exp, vector<Double_t> dexp, vector<Double_t> bkg,
              vector<Double_t> dbkg, Double_t p, Double_t s, Double_t d,
              Double_t l, Double_t Hz, Double_t Hx, Int_t lState, Int_t pState);
  bool needsARerun();
  void rerunOPModel();
  vector<Double_t> applyNormAndBackground(Double_t norm, Double_t bkg);
  /* Setters */
  void setCurrentPower(Double_t p) {currentPower = p;};
  void setCurrentStokes(Double_t s);
  void setCurrentDetune(Double_t d) {currentDetune = d;};
  void setCurrentLinewidth(Double_t l) {currentLinewidth = l;};
  void setCurrentBz(Double_t b) {currentBz = b;};
  void setCurrentBx(Double_t b) {currentBx = b;};
  void setMethod_DensityMatrix() {method = 'O';};
  void setMethod_RateEquation() {method = 'R';};
  /* Getters */
  Double_t getStokes() {return currentStokes;};

  vector<Double_t> expDat, dexpDat, expBkg, dexpBkg;   /* Data from JB, Ioana */
  vector<Double_t> mod_t, mod_e, mod_p, mod_a;         /* From my OP model */
  vector<Double_t> integrated_es;
  /*               time , fluor , pol  , align */

  /* Output files that are static */
  ClassDef(OPFitObject, 1); //Data for fitting
private:
  void updateRecents();
  Int_t getBestIndexForTime(Double_t time);
  Double_t lastPower, lastStokes, lastDetune, lastLinewidth, lastBz, lastBx;
  Double_t currentPower, currentStokes, currentDetune, currentLinewidth,
      currentBz, currentBx;
  Int_t laserState;
  Int_t polState;
  ofstream funcCalls;
  char fileName[500];
  char method;
};

ofstream funcCalls("funcCalls.txt");
