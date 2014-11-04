// Authors: Benjamin Fenker 2013
// Class that holds the data necessary to fit the optical pumping data
#include <iostream>

#include "OPFitObj.h"

#if !defined(__CINT__)
ClassImp(OPFitObject);
#endif

using std::cout;
using std::endl;

OPFitObject::OPFitObject() {
  cout << "I am here!" << endl;
}

OPFitObject::OPFitObject(vector<Double_t> exp, vector<Double_t> dexp,
                         vector<Double_t> bkg, vector<Double_t> dbkg,
                         Double_t p, Double_t s, Double_t d, Double_t l,
                         Double_t Hz, Double_t Hx, Int_t lState, Int_t pState) :
    expDat(exp), dexpDat(dexp), expBkg(bkg), dexpBkg(dbkg), mod_t(11, 0.0),
    mod_e(11, 0.0), mod_p(11, 0.0), mod_a(11, 0.0), integrated_es(11, 0.0),
    currentPower(p), currentStokes(s), currentDetune(d), currentLinewidth(l),
    currentBz(Hz), currentBx(Hx), laserState(lState), polState(pState) {
  lastPower = -1.0;
  lastStokes = s;
  lastDetune = 100.0;
  lastLinewidth = 100.0;
  lastBz = 100.0;
  lastBx = 100.0;
  funcCalls.open("funcCalls.txt", ios::trunc);

  // Get time
  time_t now;
  struct tm *tm;

  now = time(0);
  if ((tm = localtime (&now)) == NULL) {
    printf ("Error extracting time stuff\n");
    exit(1);
  }

  sprintf(fileName, "fitData%02d-%02d-%02d.tmp", tm -> tm_hour, tm -> tm_min,
          tm -> tm_sec);
}

void OPFitObject::updateRecents() {
  lastPower = currentPower;
  lastStokes = currentStokes;
  lastDetune = currentDetune;
  lastLinewidth = currentLinewidth;
  lastBz = currentBz;
  lastBx = currentBx;
}

bool OPFitObject::needsARerun() {
  return (currentPower != lastPower || currentStokes != lastStokes ||
          currentDetune != lastDetune || currentLinewidth != lastLinewidth ||
          currentBz != lastBz || currentBx != lastBx);
}

void OPFitObject::rerunOPModel() {
  char command[500];                    // Way too long!
  Double_t detune_fe = currentDetune;
  Double_t detune_ge = detune_fe - 3.6474; // From Ioana on 4/9/12
  Double_t power_fe = currentPower;
  Double_t power_ge = power_fe / 2.0;   // So says JB in an email
  
  // switch (laserState) {
  //   case 1:
  //     detune_ge = detune_fe - 1.0;
  //     break;
  //   case 2:
  //     detune_ge = detune_fe + 1.0;
  //     break;
  //   default:
  //     cout << "lState not 1 or 2.  Aborting" << endl;
  //     exit(1);
  //     break;
  // }
  sprintf(command,
          "./opticalPumping 41K 1 1740000 10.0 %g %g %18.16G %18.16G %g %18.16G %18.16G %18.16G %18.16G %s %c",
          currentBz, currentBx, detune_fe, detune_ge, currentLinewidth,
          power_fe, power_ge, currentStokes, currentStokes, fileName, method);
  // cout << command << endl;
  OPFitObject::funcCalls << command << endl;
  gSystem -> Exec(command);             // Call optical pumping code
  updateRecents();                      // Log it in baby
  ifstream opfitFile;
  opfitFile.open(fileName);
  Double_t time, es, pol, ali;
  if(opfitFile.is_open()) {
    mod_t.clear();
    mod_e.clear();
    mod_p.clear();
    mod_a.clear();
    while(!opfitFile.eof()) {
      opfitFile >> time >> es >> pol >> ali;
      mod_t.push_back(time);
      mod_e.push_back(es);
      mod_p.push_back(fabs(pol));
      mod_a.push_back(fabs(ali));
    }
  } else {
    cout << "Could not open " << fileName << " " << endl;
    exit(1);
  }
} // End rerun

vector<Double_t> OPFitObject::applyNormAndBackground(Double_t norm,
                                                     Double_t bkg) {
  vector<Double_t> res(mod_t.size(), 0.0);
  for (UInt_t i = 0; i < res.size(); i++) {
    res[i] = (norm*mod_e[i]) + bkg;
  }
  return res;
}

Int_t OPFitObject::getBestIndexForTime(Double_t time) {
  if (needsARerun()) rerunOPModel();
  Int_t trueIndex = -1;
  Double_t minTimeDiff = 1000.0;
  for (UInt_t i = 0; i < mod_t.size(); i++) {
    if (fabs(time - mod_t[i]) < minTimeDiff) {
      minTimeDiff = fabs(time - mod_t[i]);
      trueIndex = i;
    }
  }
  return trueIndex;
}

void OPFitObject::setCurrentStokes(Double_t s) {
  if (s/lastStokes < 0.0) {
    cout << "Stokes changed signs! " << endl;
    exit(1);
  }
  currentStokes = s;
}
