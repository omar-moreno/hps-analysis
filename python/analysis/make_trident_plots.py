#!/usr/bin/python

import argparse
import sys
import root_numpy as rnp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import norm

def main() : 
   
    # Use the Bayesian Methods for Hackers design
    plt.style.use('bmh')
    matplotlib.rcParams.update({'font.size': 8})

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-f", "--root_file", help="ROOT file to process.")
    args = parser.parse_args()

    if args.root_file is None : 
        print 'Failed to specify a ROOT file.'
        sys.exit(2)

    results_rec = rnp.root2array(args.root_file)

    with PdfPages("trident_plots.pdf") as pdf : 

        fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(ncols=2, nrows=2, 
                                                 figsize=(15, 15))
        ax0.hist(results_rec["vx"], bins=500, alpha=0.8, histtype="stepfilled",
                normed=True)
        ax0.set_xlim(-1, 1)

        mu = 0.013918
        sigma = 0.0844949
        ax0.set_xlim(-0.5, 0.5)
        xmin, xmax = ax0.get_xlim()
        x = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x, mu, sigma)
        ax0.plot(x, p, linewidth=2)


        ax1.hist(results_rec["vy"], bins=500, alpha=0.8, histtype="stepfilled",
                normed=True)
        #mu, sigma = norm.fit(results_rec["vy"])
        mu = -0.00363552
        sigma = 0.0209185
        ax1.set_xlim(-0.2, 0.2)
        xmin, xmax = ax1.get_xlim()
        x = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x, mu, sigma)
        ax1.plot(x, p, linewidth=2)

        ax2.hist(results_rec["vz"], bins=100, alpha=0.8, histtype="stepfilled")
        ax3.hist2d(results_rec["vx"], results_rec["vy"], bins=500)
        ax3.set_xlim(-1, 1)
        ax3.set_ylim(-0.2, 0.2)

        pdf.savefig()
        plt.close()

if __name__ == "__main__" :
    main()
