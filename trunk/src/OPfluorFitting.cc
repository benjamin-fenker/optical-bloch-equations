// Author: Benjamin Fenker

// Purpose: To write down some functions that will be useful in
// fitting OP fluorescence spectra.  
// Standard includes
#include <exception>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>

// ROOT includes
#include <TH1D.h>

// boost includes
#include <boost/algorithm/string.hpp>

namespace OPFit {
// Return value: Status (0 = success, anything else = failure)

// hist_to_fill: Initialized histogram with arbitrary binning.  It
// will be filled up by this function

// fname: Filename of OP data to fill up histogram with.  Data must
// have time in the first column and fluoresence in the last column.
// All other columsn are ignored.

// tscale: Unit scaling factor to adjust time values to match the
// histogram.  t_hist = t_file*tscale
int FillFluorHistFromFile(TH1D *hist_to_fill, std::string fname,
                          double tscale = 1.0, std::string option = "");


// Return value: whether the histograms have 1) the same number of bins and
// 2) whether the center of each bin is identical
// a, b: Histograms to check.  
bool HistsHaveSameBinning(TH1D *a, TH1D *b);

// Return value: chi2 of two histograms data, model or 0.0 if they
// have different binnigs
// data: Data histogram - errors are IGNORED
// model: model histogram - errors are USED or assume = 1 if they don't exist
// option: not used (yet?)
double CompareFluorHists(TH1D *data, TH1D *model, std::string option = "");

enum {
  success = 0,
  bad_file = 1,
  read_error = 2,
};
}

#ifndef __OPFIT__
#define __OPFIT__
int OPFit::FillFluorHistFromFile(TH1D *hist_to_fill, std::string fname,
                                 double tscale, std::string option) {
  boost::algorithm::to_lower(option);

  std::ifstream ifs(fname, std::ifstream::in);
  if (!ifs.is_open()) {
    return bad_file;
  }

  
  std::string line;
  std::vector<std::string> word;
  while (std::getline(ifs, line)) {
    while (line.c_str()[0] == ' ') line.erase(0, 1);
    //    std::cout << line << std::endl;
    boost::split(word, line, boost::is_any_of(" \t"));
    //    std::cout << word[0] << "\t" << word[word.size()-1] << std::endl;
    double time, fluor;
    try {
      time = tscale * stod(word[0]);
      fluor = stod(word[word.size()-1]);
    }
    catch (...) {
      std::cout << "Error converting to double in file " << fname << std::endl;
      return read_error;
    }
    //    std::cout << time << "\t" << fluor << std::endl;
    int binn = hist_to_fill -> FindBin(time);
    hist_to_fill ->
        SetBinContent(binn, hist_to_fill->GetBinContent(binn) + fluor);
  }
  ifs.close();
  for (int i = 1; i <= hist_to_fill -> GetNbinsX(); i++) {
    //    hist_to_fill -> SetBinError(i, sqrt(hist_to_fill->GetBinContent(i)));
    hist_to_fill -> SetBinError(i, 0.0);
    if (option == "rate") {
      hist_to_fill -> SetBinContent(i, hist_to_fill -> GetBinContent(i) / 
                                    hist_to_fill -> GetBinWidth(i));
      // hist_to_fill -> SetBinError(i, hist_to_fill -> GetBinError(i) /
      //                             hist_to_fill -> GetBinWidth(i));
    }
  }
  
  return success;
}

double OPFit::CompareFluorHists(TH1D *data, TH1D *model, std::string option) {
  int verbose = 0;
  if (verbose > 0) {
    std::cout << "Simulation\tData\tUncert.\tchi2" << std::endl;
  }

  double chi2 = 0.0;
  double t = 0.0;
  if (HistsHaveSameBinning(data, model)) {
    for (int i = 0; i <= data -> GetNbinsX(); i++) {
      double sim = model -> GetBinContent(i);
      double exp = data -> GetBinContent(i);
      double err = data -> GetBinError(i);

      t = exp - sim;
      if (err > 0.0) {
        t = t / err;
      }


      t = t * t;
      if (verbose > 0) {
        std::cout << sim << "\t" << exp << "\t" << err << "\t" << t << std::endl;
      }
      chi2 = chi2 + t;
    }
  } else {
    chi2 = 0.0;
    std::cout << "ERROR: CANNOT COMPARE HISTS " << data -> GetName() << " and ";
    std::cout << model -> GetName() << " MUST RETURN 0.0 CHI2" << std::endl;
  }
  return chi2;
}

bool OPFit::HistsHaveSameBinning(TH1D *a, TH1D *b) {
  bool same = true;
  double eps = 1.E-3;
  if (a -> GetNbinsX() != b -> GetNbinsX()) same = false;

  if (same) {
    for (int i = 1; i <= a -> GetNbinsX(); i++) {
      if (fabs(a->GetBinCenter(i) - b->GetBinCenter(i)) > eps) same = false;
    }
  }

  return same;
}

#endif                                   
