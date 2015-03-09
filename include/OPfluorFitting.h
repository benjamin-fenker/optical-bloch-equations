// Author: Benjamin Fenker
#ifndef __OPFIT__
#define __OPFIT__

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
                          double tscale = 1.0, std::string option = "",
                          double toffset = 0.);


// Return value: whether the histograms have 1) the same number of bins and
// 2) whether the center of each bin is identical
// a, b: Histograms to check.  
bool HistsHaveSameBinning(TH1D *a, TH1D *b);

// Return value: chi2 of two histograms data, model or 0.0 if they
// have different binnigs
// data: Data histogram - errors are IGNORED
// model: model histogram - errors are USED or assume = 1 if they don't exist
// option: not used (yet?)
double CompareFluorHists(TH1D *data, TH1D *model, std::string option,
                         double min, double max);

void ScaleSimulationIncludingBackground(TH1D *data, TH1D *model,
                                        double bkg_per_bin, double start,
                                        double stop);
void ScaleSimulationIncludingSNratio(TH1D *data, TH1D *model,
                                     double sn_ratio, double start,
                                     double stop);
double GetFinalPolarizationFromFile(std::string gname);
TH1D* GetResidualHistogram(TH1D *data, TH1D *model, double min, double max);

}

#endif                                   
