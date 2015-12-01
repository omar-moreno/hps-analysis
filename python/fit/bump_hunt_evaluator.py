#!/usr/bin/python

###############
#   Imports   #
###############

import argparse
import sys
import ROOT as r
from rootpy.io import root_open
import BumpHunter as bh
import numpy as np

############
#   Main   #
############

def main() : 

    # Parse all command line arguments using the argparse module
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input",  help="ROOT file containing toy histograms.")
    parser.add_argument("-n", "--number", help="Number of histograms to use in the evaluation.")
    parser.add_argument("-o", "--order",  help="Polynomial order used to model background.")
    parser.add_argument("-w", "--window", help="Window size.")
    parser.add_argument("-s", "--start",  help="Starting mass value.")
    parser.add_argument("-e", "--end",    help="End mass value.")
    args = parser.parse_args()

    if args.input is None : 
        print '[ Evaluator ]: A ROOT file needs to be specified!'
        sys.exit(2)

    # Open the ROOT file
    root_file = root_open(args.input)

    # Retrieve all the toy histograms from the file
    histos = list(root_file.objects(r.TH1))
    n_histos = len(histos)
    if args.number is not None : n_histos = int(args.number)
    print "[ Evaluator ] Using " + str(n_histos) + " histograms in evaluation of fitter."

    poly_order = 3
    if args.order is not None : poly_order = int(args.order)
    print "[ Evaluator ] Background will be fit using a " + str(poly_order) + " polynomial."

    window_size = 0.010
    if args.window is not None : window_size = float(args.window)
    print "[ Evaluator ] Using a window size " + str(window_size) 

    window_start = 0.03
    if args.start is not None : window_start = float(args.start)

    window_end = 0.04
    if args.end is not None : window_end = float(args.end)

    if (window_end - window_start) < window_size : 
        print "[ Evaluator ] Fit range is smaller than window size!"
        sys.exit(2)

    output_file_name = "eval_histos" + str(n_histos) + "_poly" + str(poly_order) + "_window" + str(window_size)
    output_file = r.TFile(output_file_name + ".root", "recreate")
    tree = r.TTree("evaluation", "evaluation")
    sig = np.zeros(1, dtype=float)
    sig_error = np.zeros(1, dtype=float)
    sig_pull = np.zeros(1, dtype=float)
    tree.Branch('signal', sig, 'signal/D')
    tree.Branch('signal_error', sig_error, 'signal_error/D')
    tree.Branch('sig_pull', sig_pull, 'sig_pull/D')

    # Prepare the fitter
    bump_hunter = bh.BumpHunter(poly_order)
    bump_hunter.set_mass_window_size(window_size)

    for histo_n in range(0, n_histos) : 
        results = bump_hunter.fit(histos[histo_n], window_start, window_end, 0.0005)
        final_params = results[0].floatParsFinal()
        sig[0] = final_params[final_params.index("nsig")].getVal()
        sig_error[0] = final_params[final_params.index("nsig")].getError()
        sig_pull[0] = sig[0]/sig_error[0]
        tree.Fill()

    output_file.Write()
    output_file.Close()

if __name__ == "__main__" : 
    main()
