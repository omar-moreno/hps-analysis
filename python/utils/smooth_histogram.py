#!/usr/bin/python

""" Script that takes a histogram as input, smooths it and applies a first order 
    interpolation to create a PDF. The resulting PDF is overlayed on the orginal
    histogram used to generate it.  The smooth histogram is saved to a ROOT file
    so it can be used to generate toys.
"""

__author__ = "Omar Moreno"
__email__  = "omoreno1@ucsc.edu"
__date__   = "February 18, 2016"

import argparse
import sys
import ROOT as r

def main() : 
   
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-f", "--root_file", help="ROOT file.")
    parser.add_argument("-n", "--name", help="Name of the histogram to smooth.")
    parser.add_argument("-i", "--iterations", 
                        help="Number of iterations to use when smoothing")
    args = parser.parse_args()

    # If a list of input files has not been specified, warn the user and exit
    # the application.
    if args.root_file is None : 
        print 'A ROOT files needs to be specified.'
        sys.exit(2)

    file = r.TFile(args.root_file)
    histogram = file.Get(args.name)

    o_file = r.TFile("smooth_histogram.root", "recreate")

    smooth_histogram = histogram.Clone()
    smooth_histogram.Smooth(int(args.iterations))

    smooth_histogram.SetAxisRange(0.07, 0.1)
    smooth_histogram.Smooth(10000, "R")
    smooth_histogram.SetAxisRange(0.0, 0.1)

    smooth_histogram.Write()

    histogram_min = histogram.GetXaxis().GetXmin()
    histogram_max = histogram.GetXaxis().GetXmax()

    x = r.RooRealVar("x", "x", 0.01, histogram_max)
    
    histogram_data = r.RooDataHist("hist_data", "hist_data", r.RooArgList(x), histogram)
    smooth_histogram_data = r.RooDataHist("smooth_hist_data", "smooth_hist_data", r.RooArgList(x), smooth_histogram)
    pdf = r.RooHistPdf("hist_pdf", "hist_pdf", r.RooArgSet(x), smooth_histogram_data, 1)
    #pdf.setUnitNorm(r.kTRUE)

    canvas = r.TCanvas("canvas", "canvas", 800, 800)

    frame = x.frame()
    histogram_data.plotOn(frame)
    pdf.plotOn(frame)
    frame.Draw()

    canvas.SetLogy()
    canvas.SaveAs("smooth_pdf.pdf")

    file.Close()
    o_file.Close()

if __name__ == "__main__" : 
    main()
    

    
