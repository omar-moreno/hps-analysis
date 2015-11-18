#!/usr/bin/python

###############
#   Imports   #
###############

import argparse
import sys
import ROOT as r

############
#   Main   #
############

def main() : 


    # Parse all command line arguments using the argparse module
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input", help="ROOT file containing the histogram.")
    parser.add_argument("-n", "--name", help="Name of histogram to generate toy from.") 
    args = parser.parse_args()

    if args.input is None : 
        print '[ Generate Toys ]: A ROOT file needs to be specified!'
        sys.exit(2)

    if args.name is None :
        print "[ Generate Toys ]: A histogram name hasn't been specified."
        sys.exit(2)

    # Open the ROOT file and retrieve the histogram of interest
    canvas = r.TCanvas("canvas", "canvas", 500, 500)
    root_file = r.TFile(args.input)
    histo = root_file.Get(args.name)
    histo.Draw()
    
    x = r.RooRealVar("x", "x", canvas.GetUxmin(), canvas.GetUxmax())
    
    #r.RooRandom.randomGenerator().SetSeed(30005);

    arg_list = r.RooArgList(x)
    histo_data = r.RooDataHist("histo_data", "histo_data", arg_list, histo)

    arg_set = r.RooArgSet(x)
    hist_pdf = r.RooHistPdf("hist_pdf", "hist_pdf", arg_set, histo_data, 0)

    generated_hist = hist_pdf.generate(arg_set, 10, r.RooFit.Extended(r.kTRUE))

    canvas.Clear()
    frame = x.frame()
    generated_hist.plotOn(frame)
    frame.Draw()

    canvas.SaveAs("test.png")

   # gauss = r.RooGaussian("gaus", "gaus", x, mean, sigma)
    #arg_set = r.RooArgSet(x)
    frame = x.frame()
    data = gauss.generate(arg_set, 1000)
    data.plotOn(frame)
    frame.Draw()

    canvas.SaveAs("test2.png")
    
    frame = x.frame()
    data = gauss.generate(arg_set, 1000)
    data.plotOn(frame)
    frame.Draw()
    
    canvas.SaveAs("test3.png")
    


if __name__ == "__main__":
    main()
