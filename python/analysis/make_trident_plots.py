#!/usr/bin/python

import argparse
import sys
import root_numpy as rnp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Ellipse
from scipy.stats import norm

def main() : 
   
    # Use the Bayesian Methods for Hackers design
    plt.style.use('bmh')
    matplotlib.rcParams.update({'font.size': 8})

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-f", "--root_file", help="ROOT file to process.")
    args = parser.parse_args()

    # If a list of input files has not been specified, warn the user and exit 
    # the app.
    if args.root_file is None : 
        print 'Failed to specify a ROOT file.'
        sys.exit(2)

    results_rec = rnp.root2array(args.root_file)

#   features = rnp.list_branches(args.root_file)
#   print '\n'.join(str(feature) for feature in features)
    
    with PdfPages("trident_plots.pdf") as pdf : 

        fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(ncols=2, nrows=2, figsize=(15, 15))

        vx_bins = np.linspace(-0.5, 0.5, 200)

        vx_arr = results_rec["vx"]
        vx_arr_cuts = vx_arr[abs(vx_arr) < 0.5]
        ax0.hist(vx_arr_cuts, vx_bins, alpha=0.8, histtype="stepfilled")
        ax0.set_xlim(-0.5, 0.5)
        
        vy_bins = np.linspace(-0.1, 0.1, 200)
        vy_arr = results_rec["vy"]
        vy_arr_cuts = vy_arr[abs(vy_arr) < 0.1]
        ax1.hist(vy_arr_cuts, vy_bins, alpha=0.8, histtype="stepfilled")
        ax1.set_xlim(-0.1, 0.1)
        
        vx_arr_cuts = vx_arr[((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1]
        ax0.hist(vx_arr_cuts, vx_bins, alpha=0.8, histtype="stepfilled")

        vy_arr_cuts = vy_arr[((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1]
        ax1.hist(vy_arr_cuts, vy_bins, alpha=0.8, histtype="stepfilled")

        ax2.hist2d(vx_arr_cuts, vy_arr_cuts, bins=300)
        ax2.set_xlim(-0.2, 0.2)
        ax2.set_ylim(-0.05, 0.05)

        #mu = 0.013918
        #sigma = 0.0844949
        #xmin, xmax = ax0.get_xlim()
        #x = np.linspace(xmin, xmax, 100)
        #p = norm.pdf(x, mu, sigma)
        #ax0.plot(x, p, linewidth=2)

        #mu = -0.00363552
        #sigma = 0.0209185
        #xmin, xmax = ax1.get_xlim()
        #x = np.linspace(xmin, xmax, 100)
        #p = norm.pdf(x, mu, sigma)
        #ax1.plot(x, p, linewidth=2)
 
        ax3.hist2d(results_rec["vx"], results_rec["vy"], bins=500)
        ax3.set_xlim(-0.5, 0.5)
        ax3.set_ylim(-0.2, 0.2)
    
        ellipse = Ellipse(xy=(0, 0), width=0.4, height=0.1, fc='None', edgecolor='r')
        ax3.add_artist(ellipse);
        
        pdf.savefig()
        plt.close()

        
        fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(ncols=2, nrows=2, figsize=(15, 15))
        p_bins = np.linspace(0, 1.5, 200)
        
        e_p_arr = results_rec["electron_p"]
        p_p_arr = results_rec["positron_p"]
        e_p_arr_cuts = e_p_arr[(e_p_arr < 1.5) & (p_p_arr < 1.5)]
        p_p_arr_cuts = p_p_arr[(e_p_arr < 1.5) & (p_p_arr < 1.5)]
        
        ax0.hist(e_p_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled")
        ax1.hist(p_p_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled")
        ax2.hist2d(e_p_arr_cuts, p_p_arr_cuts, bins=300)

        e_p_arr_cuts = e_p_arr[
               (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)
             & (e_p_arr < 1.5) 
             & (p_p_arr < 1.5)]
        p_p_arr_cuts = p_p_arr[
               (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)
             & (e_p_arr < 1.5) 
             & (p_p_arr < 1.5)]
        ax0.hist(e_p_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled")
        ax1.hist(p_p_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled")

        e_p_arr_cuts = e_p_arr[
               (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)
             & (e_p_arr < 0.85) 
             & (p_p_arr < 0.85)]
        p_p_arr_cuts = p_p_arr[
               (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)
             & (e_p_arr < 0.85) 
             & (p_p_arr < 0.85)]
        ax0.hist(e_p_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled")
        ax1.hist(p_p_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled")

        ax3.hist2d(e_p_arr_cuts, p_p_arr_cuts, bins=300)
        pdf.savefig()
        plt.close()

        bins = np.linspace(0, 1.5, 500)
        v0_p_arr = results_rec["v0_p"]
        plt.hist(v0_p_arr, bins, alpha=0.8, histtype="stepfilled")
        
        v0_p_arr_cuts = v0_p_arr[(((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled")

        v0_p_arr_cuts = v0_p_arr[
             (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)
             & (e_p_arr < 0.85) 
             & (p_p_arr < 0.85)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled")

        v0_p_arr_cuts = v0_p_arr[
             (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)
             & (e_p_arr < 0.85) 
             & (p_p_arr < 0.85)
             & (v0_p_arr > 0.85)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled")

        pdf.savefig()
        plt.close()


        bins = np.linspace(0, 0.1, 500)

        mass_arr = results_rec["invariant_mass"]
        mass_arr_cuts = mass_arr[mass_arr < 0.1]
        plt.hist(mass_arr, bins, alpha=0.8, histtype="stepfilled")
        plt.xlabel("Invariant Mass M(e^-e^+) (GeV)")
        plt.ylabel("Events/.2")
        plt.xlim(0, 0.1)
       
        mass_arr_cuts = mass_arr[((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1]
        plt.hist(mass_arr_cuts, bins, alpha=0.8, histtype="stepfilled") 

        mass_arr_cuts = mass_arr[
               (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)
             & (e_p_arr < 0.85) 
             & (p_p_arr < 0.85)]
        plt.hist(mass_arr_cuts, bins, alpha=0.8, histtype="stepfilled") 

        mass_arr_cuts = mass_arr[
               (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)
             & (e_p_arr < 0.85) 
             & (p_p_arr < 0.85)
             & (v0_p_arr > 0.85)]
        plt.hist(mass_arr_cuts, bins, alpha=0.8, histtype="stepfilled") 
        #mass_arr = mass_arr[v0_p > 0.8]

        pdf.savefig()
        plt.close()


if __name__ == "__main__" :
    main()
