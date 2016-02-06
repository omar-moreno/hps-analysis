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

    # Parse all command line arguments 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--ideal",  action='store_true',  
                        help="Use an ideal distribution to generate toys") 
    parser.add_argument("-f", "--file",
                        help="ROOT file containing the histogram.")
    parser.add_argument("-n", "--name", 
                        help="Name of histogram to generate toys from.") 
    parser.add_argument("-c", "--toy_count",
                        help="Number of toy histograms to generate.")
    parser.add_argument("-e", "--events",
                        help="Number of events per histogram.")
    parser.add_argument("-o", "--output",
                        help="ROOT output file.")
    parser.add_argument("-b", "--bins",
                        help="Number of bins the histgram should have.")
    args = parser.parse_args()

    # If not using an ideal distribution, check that a file and histogram name
    # were specified.
    root_file = None
    histo = None
    if not args.ideal : 

        if args.file is None :
            print '[ Generate Toys ]: Please specify a ROOT file to process.'
            sys.exit(2)
        elif args.name is None :
            print "[ Generate Toys ]: A histogram name wasn't specified."
            sys.exit(2)
        else :
            root_file = r.TFile(args.file)
            histo = root_file.Get(args.name)

    else : print "[ Generate Toys ]: Using ideal distribution."

    # Set the default number of toy histograms to generate
    toy_count = 1
    if args.toy_count is not None: toy_count = int(args.toy_count)
    print '[ Generate Toys ]: Creating ' + str(toy_count) + ' histograms.'

    events = 1000
    if args.events is not None : events = int(args.events)
    print '[ Generate Toys ]: Generating ' + str(events) + ' per histogram.'

    #
    x = None
    if args.ideal : 
        x = r.RooRealVar("x", "x", 171, 264)
    else : 
        x = r.RooRealVar("x", "x", 0, 0.1)
    if args.bins is not None : x.setBins(int(args.bins))

    pdf = None
    if args.ideal : 
        pdf = r.RooGenericPdf("apex_pdf", "apex_pdf", "pow(170-x, 2)*pow(265-x,2)/pow(x,4)", r.RooArgList(x))
    else : 
        histo_data = r.RooDataHist("histo_data", "histo_data", r.RooArgList(x), histo)
        pdf = r.RooHistPdf("hist_pdf", "hist_pdf", r.RooArgSet(x), histo_data, 7)


    canvas = r.TCanvas("canvas", "canvas", 800, 800)
    frame = x.frame()
    pdf.plotOn(frame)
    frame.Draw()
    canvas.SaveAs("pdf.pdf")


    output_file_name = "output.root"
    if args.output is not None: output_file_name = args.output
    output_file = r.TFile(output_file_name, "RECREATE")

    for hist_count in range(0, toy_count) : 
        generated_hist = pdf.generateBinned(r.RooArgSet(x), events, r.RooFit.Extended(r.kTRUE))
        generated_hist.createHistogram("toy_" + str(hist_count), x).Write()

    output_file.Close()

if __name__ == "__main__":
    main()
