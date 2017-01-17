#!/usr/bin/python

import argparse
import sys
import root_numpy as rnp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import norm

def plot_sig(pdf, mass, sig, sig_err, bins, **params):
        
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3, figsize=(30, 10))
    fig.suptitle("A' mass: " + str(mass) + " GeV")

    label=None
    if 'label' in params:
        label=params['label']
        ax.legend()
        
    ax0.hist(sig, bins, histtype='step', lw=1.5, label=label)
    ax0.set_xlabel('Signal Yield')
    ax0.set_ylabel('AU')

    ax1.hist(sig_err, bins, histtype='step', lw=1.5, label=label)
    ax1.set_xlabel('Signal Yield Error')
    ax1.set_ylabel('AU')

    pull = sig/sig_err
    ax2.hist(pull, bins, histtype='step', lw=1.5, label=label)
    ax2.set_xlabel('Pull')
    ax2.set_ylabel('AU')

    pdf.savefig(bbox_inches='tight')
    plt.close()

def plot_bkg(pdf, mass, bkg, bkg_err, bins, **params):
        
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3, figsize=(30, 10))
    fig.suptitle("A' mass: " + str(mass) + " GeV")

    label=None
    if 'label' in params:
        label=params['label']
        ax.legend()
        
    ax0.hist(bkg, bins, histtype='step', lw=1.5, label=label)
    ax0.set_xlabel('Background Yield')
    ax0.set_ylabel('AU')

    ax1.hist(bkg_err, bins, histtype='step', lw=1.5, label=label)
    ax1.set_xlabel('Background Yield Error')
    ax1.set_ylabel('AU')

    pull = bkg/bkg_err
    ax2.hist(pull, bins, histtype='step', lw=1.5, label=label)
    ax2.set_xlabel('Pull')
    ax2.set_ylabel('AU')

    pdf.savefig(bbox_inches='tight')
    plt.close()

def main(): 
   
    # Use the Bayesian Methods for Hackers design
    plt.style.use('bmh')
    matplotlib.rcParams.update({'font.size': 8})

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-l", "--file_list", help="ROOT file")
    args = parser.parse_args()

    # If a list of files has not been specified, warn the user and exit
    # the application.
    if not args.file_list: 
        print 'A list of ROOT files needs to be specified'
        sys.exit(2)

    # Open the file containing the list of files to process
    root_file_list = None
    try:
        root_file_list = open(args.file_list, 'r')
    except IOError: 
        print 'Unable to open file %s' % args.file_list
        sys.exit(2)

    root_files = []
    for line in root_file_list: 
        root_files.append(line.strip())

    rec = rnp.root2array(root_files, 'results')

    make_plots(rec)

def make_plots(rec):

    # Use the 'Bayesian Methods for Hackers' style
    plt.style.use('bmh')
    matplotlib.rcParams.update({'font.size': 12})
    matplotlib.rcParams['axes.facecolor'] = 'white'
    matplotlib.rcParams['legend.numpoints'] = 1        

    # Set the default font size for plots and legends    
    matplotlib.rcParams.update({'font.size': 12, 'legend.fontsize': 9})

    ap_mass = rec['ap_mass']
    bkg_yield = rec['bkg_yield']
    bkg_yield_err = rec['bkg_yield_err']
    sig_yield = rec['sig_yield']
    sig_yield_err = rec['sig_yield_err']

    bkg_mean = []
    bkg_pull = []
    sig_mean = []
    sig_pull = []

    # Create a unique list of A' masses
    ap_masses = np.unique(ap_mass)

    with PdfPages("results.pdf") as pdf :

        for mass in ap_masses :
            
            print "Processing A' mass: " + str(mass)

            bkg = bkg_yield[ap_mass == mass]
            bkg_err = bkg_yield_err[ap_mass == mass]
            plot_bkg(pdf, mass, bkg, bkg_err, 15)
            bkg_mean.append(np.mean(bkg))
            bkg_pull.append(np.mean((bkg - np.mean(bkg))/bkg_err))

            sig = sig_yield[ap_mass == mass]
            sig_err = sig_yield_err[ap_mass == mass]
            plot_sig(pdf, mass, sig, sig_err, 15)
            sig_mean.append(np.mean(sig))
            sig_pull.append(np.mean(sig/sig_err))
       
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 10))
        ax.plot(ap_masses, bkg_mean, 'o-')
        plt.xlabel("$A'$ mass hypothesis")
        plt.ylabel("Background Yield")
        pdf.savefig()
        plt.close()

        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 10))
        ax.plot(ap_masses, bkg_pull, 'o-')
        plt.xlabel("$A'$ mass hypothesis")
        plt.ylabel("Background Yield Pull")
        plt.ylim([-2, 2])
        pdf.savefig()
        plt.close()

        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 10))
        ax.plot(ap_masses, sig_mean, 'o-')
        plt.xlabel("$A'$ mass hypothesis")
        plt.ylabel("Signal Yield")
        pdf.savefig()
        plt.close()

        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 10))
        ax.plot(ap_masses, sig_pull, 'o-')
        plt.xlabel("$A'$ mass hypothesis")
        plt.ylabel("Signal Yield Pull")
        plt.ylim([-2, 2])
        pdf.savefig()
        plt.close()


if __name__ == "__main__" : 
    main()

