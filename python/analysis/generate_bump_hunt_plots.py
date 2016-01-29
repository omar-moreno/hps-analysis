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
    parser.add_argument("-l", "--input_list", help="ROOT file")
    parser.add_argument("-b", "--bkg_only", help="Process a file generated with a bkg only fit.")
    args = parser.parse_args()

    # If a list of input files has not been specified, warn the user and exit
    # the application.
    if args.input_list is None : 
        print 'A list of ROOT files needs to be specified'
        sys.exit(2)

    # Open the file containing the list of files to process
    try:
        file_list = open(args.input_list, 'r')
    except IOError: 
        print "Unable to open file " + str(args.input_list)
        sys.exit(2)

    # Create a map between the name of a variable in the tuple and its position
    # along the array.
    index_map = {
            "ap_mass"           : 0, 
            "sig_yield"         : 1, 
            "sig_yield_error"   : 2, 
            "bkg_yield"         : 3, 
            "bkg_yield_error"   : 4, 
            "nll"               : 5, 
            "invalid_nll"       : 6, 
            "minuit_status"     : 7, 
            "edm"               : 8
    }

    variable_map = { 
        "bkg_mean"  : [],
        "bkg_sigma" : [],
        "bkg_pull"  : []

    }

    name_list = []

    results_list = []
    yield_pull_array_list = []

    for line in file_list : 
        
        print "Processing file: " + str(line.strip())

        results_rec = rnp.root2array(line.strip())

        results_list.append(rnp.rec2array(results_rec))
        
        variable_map["bkg_mean"].append([])
        variable_map["bkg_sigma"].append([])
        variable_map["bkg_pull"].append([])

        split_string = line.split('_', 2)
        #print split_string
        name_list.append("order " + str(split_string[0][5]))

        yield_pull_array_list.append([])


    # Create a unique list of A' masses
    ap_masses = np.unique(results_list[0][:,0])

    with PdfPages("results.pdf") as pdf :

        for ap_mass in ap_masses :
            
            print "Processing A' mass: " + str(ap_mass)
            
            fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, figsize=(15, 5))
            
            index = 0
            for results in results_list : 
            
                bkg_yield_a = results[:,index_map["bkg_yield"]][results[:,index_map["ap_mass"]] == ap_mass]
                bkg_yield_err_a = results[:,index_map["bkg_yield_error"]][results[:,index_map["ap_mass"]] == ap_mass]
            
                ax0.hist(bkg_yield_a, bins=50, alpha=0.8, histtype="stepfilled", normed=True, label=name_list[index])
                ax0.set_xlabel('Background Yield', fontsize=8)
                ax0.set_ylabel('AU', fontsize=8)
                ax0.set_title("$A'$ mass hypothesis = " + str(ap_mass), fontsize=8)
                ax0.legend()

                mu, sigma = norm.fit(bkg_yield_a) 
                variable_map["bkg_mean"][index].append(mu)
                variable_map["bkg_sigma"][index].append(sigma)

                xmin, xmax = ax0.get_xlim()
                x = np.linspace(xmin, xmax, 100)
                p = norm.pdf(x, mu, sigma)
                ax0.plot(x, p, linewidth=2)

                ax1.hist(bkg_yield_err_a, bins=50, alpha=0.8, histtype="stepfilled", normed=True, label=name_list[index])
                ax1.set_xlabel('Background Yield Fit Error', fontsize=8)
                ax1.set_ylabel('AU', fontsize=8)
                ax1.set_title("$A'$ mass hypothesis = " + str(ap_mass), fontsize=8)
                ax1.legend()

                ax2.hist((bkg_yield_a - mu)/bkg_yield_err_a, bins=50, alpha=0.8, histtype="stepfilled", normed=True, label=name_list[index])
                ax2.set_xlabel('Background Yield Pull', fontsize=8)
                ax2.set_ylabel('AU', fontsize=8)
                ax2.set_title("$A'$ mass hypothesis = " + str(ap_mass),fontsize=8)
                ax2.legend()

                mu, sigma = norm.fit((bkg_yield_a - mu)/bkg_yield_err_a)
                variable_map["bkg_pull"][index].append(mu)

                xmin, xmax = ax2.get_xlim()
                x = np.linspace(xmin, xmax, 100)
                p = norm.pdf(x, mu, sigma)
                ax2.plot(x, p, linewidth=2)

                index += 1

            pdf.savefig()
            plt.close()
        
        for index in range(0, len(variable_map["bkg_mean"])) : 

            bkg_yield_mean_a = np.array(variable_map["bkg_mean"][index])

            plt.plot(ap_masses, variable_map["bkg_mean"][index], 'o-', label=name_list[index])
            plt.xlabel("$A'$ mass hypothesis")
            plt.ylabel("Background Yield")
            plt.legend()
            #ax0.fill_between(ap_masses, yield_mean_array - yield_sigma_array, yield_mean_array + yield_sigma_array)
            
        pdf.savefig()
        plt.close()

        for index in range(0, len(variable_map["bkg_mean"])) : 
        
            plt.plot(ap_masses, variable_map["bkg_pull"][index], 'o-', label=name_list[index])
            plt.xlabel("$A'$ mass hypothesis")
            plt.ylabel("Background Yield Pull")
            plt.ylim([-0.5, 0.5])
            plt.legend()

        pdf.savefig()
        plt.close()

if __name__ == "__main__" : 
    main()

