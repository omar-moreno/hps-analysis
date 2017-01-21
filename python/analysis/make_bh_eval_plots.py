#!/usr/bin/python

import argparse
import sys
import root_numpy as rnp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import norm

def plot_sig(axes, sig, sig_err, bins, **params):

    ax = axes.flatten()

    label=''
    if 'label' in params:
        label = params['label']

    ax[0].hist(sig, bins, histtype='step', lw=1.5, label=label)
    ax[0].set_xlabel('Signal Yield')
    ax[0].set_ylabel('AU')
    ax[0].legend()

    ax[1].hist(sig_err, bins, histtype='step', lw=1.5, label=label)
    ax[1].set_xlabel('Signal Yield Error')
    ax[1].set_ylabel('AU')
    ax[1].legend()

    pull = sig/sig_err
    ax[2].hist(pull, bins, histtype='step', lw=1.5, label=label)
    ax[2].set_xlabel('Pull')
    ax[2].set_ylabel('AU')
    ax[2].legend()

def plot_bkg(axes, bkg, bkg_err, bins, **params):

    ax = axes.flatten()

    label=''
    if 'label' in params:
        label=params['label']

    ax[0].hist(bkg, bins, histtype='step', lw=1.5, label=label)
    ax[0].set_xlabel('Background Yield')
    ax[0].set_ylabel('AU')
    ax[0].legend()

    ax[1].hist(bkg_err, bins, histtype='step', lw=1.5, label=label)
    ax[1].set_xlabel('Background Yield Error')
    ax[1].set_ylabel('AU')
    ax[1].legend()

    pull = bkg/bkg_err
    ax[2].hist(pull, bins, histtype='step', lw=1.5, label=label)
    ax[2].set_xlabel('Pull')
    ax[2].set_ylabel('AU')
    ax[2].legend()

def plot_ul(axes, ul, bins, label):
    
    axes.hist(ul, bins, histtype='step', lw=1.5, label=label)
    axes.set_xlabel('Signal Upper Limit')
    axes.set_ylabel('AU')
    axes.legend()

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

    ap_mass    = rec['ap_mass']
    poly_order = rec['poly_order']

    bkg_yield     = rec['bkg_yield']
    bkg_yield_err = rec['bkg_yield_err']

    sig_yield     = rec['sig_yield']
    sig_yield_err = rec['sig_yield_err']

    upper_limit = rec['upper_limit']

    # Create a unique list of A' masses
    ap_masses = np.unique(ap_mass)

    # Create a unique list of Polynomials
    polys = np.unique(poly_order)

    bkg_mean = {}
    bkg_pull = {}
    sig_mean = {}
    sig_pull = {}
    ul_mean = {}

    bkg_figures = {}
    bkg_axes = {}
    sig_figures = {}
    sig_axes = {}
    ul_figures = {}
    ul_axes = {}
    for mass in ap_masses:
        bkg_figures[mass], bkg_axes[mass] = plt.subplots(nrows=1, ncols=3, figsize=(30, 10))
        bkg_figures[mass].suptitle("A' mass: %s GeV" %mass)
        
        sig_figures[mass], sig_axes[mass] = plt.subplots(nrows=1, ncols=3, figsize=(30, 10))
        sig_figures[mass].suptitle("A' mass: %s GeV" % mass)

        ul_figures[mass], ul_axes[mass] = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
        ul_figures[mass].suptitle("A' mass: %s GeV" % mass)  

    poly_fig = {}
    poly_ax = {}
    for fig_n in xrange(0, 5):
        poly_fig[fig_n], poly_ax[fig_n] = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

    with PdfPages("results.pdf") as pdf :

        for poly in polys:
            
            bkg_mean[poly] = []
            bkg_pull[poly] = []
            sig_mean[poly] = []
            sig_pull[poly] = []
            ul_mean[poly]  = []

            for mass in ap_masses :

                selection = (ap_mass == mass) & (poly_order == poly)

                print "Processing A' mass: " + str(mass)

                bkg = bkg_yield[selection]
                bkg_err = bkg_yield_err[selection]
                plot_bkg(bkg_axes[mass], bkg, bkg_err, 15, label=("poly order: %s" % poly))
                bkg_mean[poly].append(np.mean(bkg))
                bkg_pull[poly].append(np.mean((bkg - np.mean(bkg))/bkg_err))

                sig = sig_yield[selection]
                sig_err = sig_yield_err[selection]
                plot_sig(sig_axes[mass], sig, sig_err, 15, label=("poly order: %s" % poly))
                sig_mean[poly].append(np.mean(sig))
                sig_pull[poly].append(np.mean(sig/sig_err))

                plot_ul(ul_axes[mass], upper_limit[selection], 15, ("poly order: %s" % poly))
                ul_mean[poly].append(np.mean(upper_limit[selection]))

            poly_ax[0].plot(ap_masses, bkg_mean[poly], 'o-', label=("poly order: %s" % poly))
            poly_ax[0].set_xlabel("$A'$ mass hypothesis")
            poly_ax[0].set_ylabel("Background Yield")

            poly_ax[1].plot(ap_masses, bkg_pull[poly], 'o-', label=("poly order: %s" % poly))
            poly_ax[1].set_xlabel("$A'$ mass hypothesis")
            poly_ax[1].set_ylabel("Background Yield Pull")
            poly_ax[1].set_ylim([-2, 2])

            poly_ax[2].plot(ap_masses, sig_mean[poly], 'o-', label=("poly order: %s" % poly))
            poly_ax[2].set_xlabel("$A'$ mass hypothesis")
            poly_ax[2].set_ylabel("Signal Yield")

            poly_ax[3].plot(ap_masses, sig_pull[poly], 'o-', label=("poly order: %s" % poly))
            poly_ax[3].set_xlabel("$A'$ mass hypothesis")
            poly_ax[3].set_ylabel("Signal Yield Pull")
            poly_ax[3].set_ylim([-2, 2])

            poly_ax[4].plot(ap_masses, ul_mean[poly], 'o-', label=("poly order: %s" % poly))
            poly_ax[4].set_xlabel("$A'$ mass hypothesis")
            poly_ax[4].set_ylabel("Signal Upper Limit")
        
        for mass in ap_masses:
            pdf.savefig(bkg_figures[mass], bbox_inches='tight')
            pdf.savefig(sig_figures[mass], bbox_inches='tight')
            pdf.savefig(ul_figures[mass], bbox_inches='tight')

        for fig_n in xrange(0, 5):
            pdf.savefig(poly_fig[fig_n], bbox_inches='tight')

if __name__ == "__main__" : 
    main()

