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

    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input", help="ROOT file")
    args = parser.parse_args()

    if args.input is None : 
        print 'A ROOT file needs to be specified'
        sys.exit(2)

    results_rec = rnp.root2array(args.input)
    results_array = rnp.rec2array(results_rec)

    ap_masses = np.unique(results_array[:,0])
    yield_mean_list = []
    yield_sigma_list = []
    yield_pull_list = []

    with PdfPages("results.pdf") as pdf :

        for ap_mass in ap_masses : 
            print "Processing A' mass: " + str(ap_mass)
            
            yield_array = results_array[:,1][results_array[:,0] == ap_mass]
            yield_error_array = results_array[:,2][results_array[:,0] == ap_mass]
            
            fig, (ax0, ax1, ax2) = plt.subplots(ncols=3)

            ax0.hist(yield_array, bins=100, alpha=0.8, histtype="stepfilled", normed=True)
            ax0.set_xlabel('Background yield', fontsize=8)
            ax0.set_ylabel('AU', fontsize=8)
            ax0.set_title("$A'$ mass hypothesis = " + str(ap_mass), fontsize=8)

            mu, std = norm.fit(yield_array)
            yield_mean_list.append(mu)
            yield_sigma_list.append(std)
            xmin, xmax = ax0.get_xlim()
            x = np.linspace(xmin, xmax, 100)
            p = norm.pdf(x, mu, std)
            ax0.plot(x, p, linewidth=2)

            ax1.hist(yield_error_array, bins=100, alpha=0.8, histtype="stepfilled", normed=True)
            ax1.set_xlabel('Background yield fit error', fontsize=8)
            ax1.set_ylabel('AU', fontsize=8)
            ax1.set_title("$A'$ mass hypothesis = " + str(ap_mass), fontsize=8)

            ax2.hist(yield_array/yield_error_array, bins=100, alpha=0.8, histtype="stepfilled", normed=True)
            ax2.set_xlabel('Pull', fontsize=8)
            ax2.set_ylabel('AU', fontsize=8)
            ax2.set_title("$A'$ mass hypothesis = " + str(ap_mass),fontsize=8)

            mu, std = norm.fit(yield_array/yield_error_array)
            yield_pull_list.append(mu)
            xmin, xmax = ax2.get_xlim()
            x = np.linspace(xmin, xmax, 100)
            p = norm.pdf(x, mu, std)
            ax2.plot(x, p, linewidth=2)

            pdf.savefig()
            plt.close()
        
        yield_mean_array = np.array(yield_mean_list)
        yield_sigma_array = np.array(yield_sigma_list)
        yield_pull_array = np.array(yield_pull_list)

        fig, (ax0, ax1) = plt.subplots(ncols=2)
        ax0.plot(ap_masses, yield_mean_array, 'o-')
        ax0.set_xlabel("$A'$ mass hypothesis")
        ax0.set_ylabel("Background yield")
        ax0.fill_between(ap_masses, yield_mean_array - yield_sigma_array, yield_mean_array + yield_sigma_array)

        ax1.plot(ap_masses, yield_pull_array, 'o-')
        ax1.set_xlabel("$A'$ mass hypothesis")
        ax1.set_ylabel("Background yield pull")

        pdf.savefig()
        plt.close()

if __name__ == "__main__" : 
    main()

