// Author: Benjamin Fenker
#include <iostream>
#include <string>
#include <vector>

#include <exception>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>

// ROOT includes
#include <TH1D.h>

// boost includes
#include <boost/algorithm/string.hpp>

#include "OPfluorFitting.h"
enum {
  success = 0,
  bad_file = 1,
  read_error = 2,
};


int OPFit::FillFluorHistFromFile(TH1D *hist_to_fill, std::string fname,
                                 double tscale, std::string option,
                                 double toffset) {
  boost::algorithm::to_lower(option);

  std::ifstream ifs(fname, std::ifstream::in);
  if (!ifs.is_open()) {
    return bad_file;
  }

  
  std::string line;
  std::vector<std::string> word;
  std::vector<int> times_filled(hist_to_fill -> GetNbinsX() + 1, 0);
  for (int i = 1; i <= hist_to_fill -> GetNbinsX(); i++) 
    hist_to_fill -> SetBinContent(i, 0.0);

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
    // Adjust time
    time = time + toffset;
    //    std::cout << time << "\t" << fluor << std::endl;
    int binn = hist_to_fill -> FindBin(time);
    hist_to_fill ->
        SetBinContent(binn, hist_to_fill->GetBinContent(binn) + fluor);
    times_filled[binn]++;
  }
  ifs.close();
  for (int i = 1; i <= hist_to_fill -> GetNbinsX(); i++) {
    if (times_filled[i] > 0) {
      hist_to_fill ->
          SetBinContent(i, hist_to_fill -> GetBinContent(i) / times_filled[i]);
    } else {
      hist_to_fill -> SetBinContent(i, 0.0);
    }
    //    hist_to_fill -> SetBinContent(i, hist_to_fill -> GetBinContent(i));
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

TH1D* OPFit::GetResidualHistogram(TH1D *data, TH1D *model,
                                  double min, double max) {
  if (!HistsHaveSameBinning(data, model)) {
    std::cout << "ERROR: CANNOT COMPARE HISTS " << data -> GetName() << " and ";
    std::cout << model -> GetName() << " MUST RETURN 0.0 CHI2" << std::endl;
    return NULL;
  }
  TH1D *residuals = new TH1D("residuals", "residuals", data -> GetNbinsX(),
                             data -> GetBinLowEdge(1),
                             data -> GetBinLowEdge(data -> GetNbinsX()) +
                             data -> GetBinWidth(data -> GetNbinsX()));
  for (int i = 0; i <= data -> GetNbinsX(); i++) {
      if (data -> GetBinLowEdge(i) < min) continue;
      if (data -> GetBinLowEdge(i) + data -> GetBinWidth(i) > max) continue;

      double sim = model -> GetBinContent(i);
      double exp = data -> GetBinContent(i);
      double r = 0.0;
      if (sim > 0) {
        r = (exp - sim) / sqrt(sim);
      }
      residuals -> SetBinContent(i, r);
      residuals -> SetBinError(i, 0.0);
  }
  return residuals;
}

double OPFit::CompareFluorHists(TH1D *data, TH1D *model, std::string option,
                                double min, double max) {
  int verbose = 0;
  if (verbose > 0) {
    std::cout << "Center\tSimulation\tData\tUncert.\tchi2" << std::endl;
  }
  boost::algorithm::to_lower(option);
  if (!HistsHaveSameBinning(data, model)) {
    std::cout << "ERROR: CANNOT COMPARE HISTS " << data -> GetName() << " and ";
    std::cout << model -> GetName() << " MUST RETURN 0.0 CHI2" << std::endl;
    return 0.0;
  }
  if (option == "logl") {
    //    std::cout << "Fitting logl" << std::endl;
    double logl = 0.0;
    double t = 0.0;
    for (int i = 0; i <= data -> GetNbinsX(); i++) {
      if (data -> GetBinLowEdge(i) < min) continue;
      if (data -> GetBinLowEdge(i) + data -> GetBinWidth(i) > max) continue;
      double sim = model -> GetBinContent(i);
      double exp = data -> GetBinContent(i);

      if (sim > 0) {
        if (exp > 0) {
          t = (sim - exp + exp*log(exp / sim));
        } else {                        // exp <= 0
          t = sim;
        }
      } else {                          //  sim < 0
        t = 0.0;
      }
      if (verbose > 0) {
        std::cout << data -> GetBinCenter(i) << "\t" << sim << "\t"
                  << exp << "\t" << t << std::endl;
      }

      // Multiply by two to get error definitions right
      t = t * 2.0;

      logl = logl + t;
    }
    int j;
    if (verbose > 0) {
      std::cout << "Enter any number to continue...";
      std::cin >> j;
    }
    return logl;
  } else {
    double chi2 = 0.0;
    double t = 0.0;
    for (int i = 0; i <= data -> GetNbinsX(); i++) {
      
      if (data -> GetBinLowEdge(i) < min) continue;
      if (data -> GetBinLowEdge(i) + data -> GetBinWidth(i) > max) continue;
      double sim = model -> GetBinContent(i);
      double exp = data -> GetBinContent(i);
      double err = data -> GetBinError(i);

      t = exp - sim;
      if (err > 0.0) {
        t = t / (err*err);              // squared right? of course!
      }
      
      t = t * t;
      if (verbose > 0) {
        std::cout << data -> GetBinCenter(i) << "\t" << sim << "\t"
                  << exp << "\t" << err << "\t" << t << std::endl;
      }
      chi2 = chi2 + t;
    }
    return chi2;
  }
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

void OPFit::ScaleSimulationIncludingBackground(TH1D *data, TH1D *model,
                                               double bkg_per_bin, double start,
                                               double stop) {
  double scale_num = (data -> Integral(data -> FindBin(start),
                                       data -> FindBin(stop)));
  int nbins = data -> FindBin(stop) - data -> FindBin(start);
  // cout << "In the data there are " << scale_num << " counts in "
  //      << nbins << " bins" << endl;

  scale_num = scale_num - (bkg_per_bin * nbins);
  double scale_den = model -> Integral(model -> FindBin(start),
                                       model -> FindBin(stop));
      
  double scale = scale_num / scale_den;
  //  scale = 2633.;
  //  std::cout << " Scale factor = " << scale << " ---> ";
  model -> Scale(scale);
  for (int i = 1; i < model -> GetNbinsX(); i++) {
    model -> SetBinContent(i, model -> GetBinContent(i) + bkg_per_bin);
  }

}

void OPFit::ScaleSimulationIncludingSNratio(TH1D *data, TH1D *model,
                                            double sn_ratio, double start,
                                            double stop) {
  int nbins = data -> FindBin(stop) - data -> FindBin(start);
  double ntot = (data -> Integral(data -> FindBin(start),
                                  data -> FindBin(stop)));
  double bkg_per_bin = (ntot / nbins) / (sn_ratio + 1.0);
  // std::cout << "Ntot: " << ntot << " nbins: " << nbins << " sn_ratio: "
  //           << sn_ratio << " bkg_per_bin: " << bkg_per_bin << std::endl;

  ScaleSimulationIncludingBackground(data, model, bkg_per_bin, start, stop);
}

double OPFit::GetFinalPolarizationFromFile(std::string fname) {
  std::ifstream ifs(fname, std::ifstream::in);
  if (!ifs.is_open()) {
    std::cout << "Could not open " << fname << std::endl;
    return 0.0;
  }
  std::string line;
  std::vector<std::string> word;
  while (std::getline(ifs, line)) {
    boost::split(word, line, boost::is_any_of(" \t"));
  }
  double polarization = stod(word[word.size() - 3]);
  std::cout << "Polarization = " << polarization << std::endl;
  return polarization;
  
}
