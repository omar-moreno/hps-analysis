#!/usr/bin/python

import argparse
import sys
import root_numpy as rnp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import norm

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
    matplotlib.rcParams.update({'font.size': 16, 'legend.fontsize': 9})

    # Enable the use of LaTeX in titles
    plt.rc('text', usetex=True)
  
    ap_mass = rec['ap_mass']
    sig_yield = rec['sig_yield']
    sig_yield_err = rec['sig_yield_err']

    # Create a unique list of A' masses
    ap_masses = np.unique(ap_mass)

    with PdfPages("results.pdf") as pdf :

        for mass in ap_masses :
            
            print "Processing A' mass: " + str(mass)
            
            fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, figsize=(15, 5))
            ax0.hist(sig_yield[ap_mass == mass], bins=50, histtype="step", normed=True)
            ax0.set_xlabel('signal yield')
            ax0.set_ylabel('au')
            ax0.legend()
        
            mu, sigma = norm.fit(sig_yield[ap_mass == mass]) 
            xmin, xmax = ax0.get_xlim()
            x = np.linspace(xmin, xmax, 100)
            p = norm.pdf(x, mu, sigma)
            ax0.plot(x, p, linewidth=2)

            ax1.hist(sig_yield_err[ap_mass == mass], bins=50, histtype="step", normed=True)
            ax1.set_xlabel('Signal Yield Fit Error')
            ax1.set_ylabel('AU')
            ax1.legend()

            pull = sig_yield[ap_mass == mass]/sig_yield_err[ap_mass == mass]
            ax2.hist(pull, bins=50, histtype="step", normed=True)
            ax2.set_xlabel('Signal Yield Pull')
            ax2.set_ylabel('AU')
           
            pdf.savefig()
            plt.close()

'''
            fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, figsize=(15, 5))
            
            index = 0
            for results in results_list : 
            
                bkg_yield_a = results["bkg_yield"]
                bkg_yield_a = bkg_yield_a[results["ap_mass"] == ap_mass]
                bkg_yield_err_a = results["bkg_yield_err"]
                bkg_yield_err_a = bkg_yield_err_a[results["ap_mass"] == ap_mass]
            
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
                '''

'''
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

            '''
'''
            pdf.savefig()
            plt.close()
            
            if not args.bkg_only : 

                fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, figsize=(15, 5))
            
                index = 0
                for results in results_list : 
            
                    sig_yield_a = results["sig_yield"][results["ap_mass"] == ap_mass]
                    sig_yield_err_a = results["sig_yield_err"][results["ap_mass"] == ap_mass]
            
                    ax0.hist(sig_yield_a, bins=50, alpha=0.8, histtype="stepfilled", normed=true, label=name_list[index])
                    ax0.set_xlabel('signal yield', fontsize=8)
                    ax0.set_ylabel('au', fontsize=8)
                    ax0.set_title("$a'$ mass hypothesis = " + str(ap_mass), fontsize=8)
                    ax0.legend()

                    mu, sigma = norm.fit(sig_yield_a) 
                    variable_map["sig_mean"][index].append(mu)
                    variable_map["sig_sigma"][index].append(sigma)

                    xmin, xmax = ax0.get_xlim()
                    x = np.linspace(xmin, xmax, 100)
                    p = norm.pdf(x, mu, sigma)
                    ax0.plot(x, p, linewidth=2)

                    ax1.hist(sig_yield_err_a, bins=50, alpha=0.8, histtype="stepfilled", normed=True, label=name_list[index])
                    ax1.set_xlabel('Signal Yield Fit Error', fontsize=8)
                    ax1.set_ylabel('AU', fontsize=8)
                    ax1.set_title("$A'$ mass hypothesis = " + str(ap_mass), fontsize=8)
                    ax1.legend()

                    ax2.hist(sig_yield_a/sig_yield_err_a, bins=50, alpha=0.8, histtype="stepfilled", normed=True, label=name_list[index])
                    ax2.set_xlabel('Signal Yield Pull', fontsize=8)
                    ax2.set_ylabel('AU', fontsize=8)
                    ax2.set_title("$A'$ mass hypothesis = " + str(ap_mass),fontsize=8)
                    ax2.legend()

                    mu, sigma = norm.fit(sig_yield_a/sig_yield_err_a)
                    variable_map["sig_pull"][index].append(mu)

                    xmin, xmax = ax2.get_xlim()
                    x = np.linspace(xmin, xmax, 100)
                    p = norm.pdf(x, mu, sigma)
                    ax2.plot(x, p, linewidth=2)

                    index += 1

                pdf.savefig()
                plt.close()

        for index in range(0, len(variable_map["bkg_mean"])) : 

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

        if not args.bkg_only : 
        
            for index in range(0, len(variable_map["sig_mean"])) : 

                plt.plot(ap_masses, variable_map["sig_mean"][index], 'o-', label=name_list[index])
                plt.xlabel("$A'$ mass hypothesis")
                plt.ylabel("Signal Yield")
                plt.legend()
                #plt.fill_between(ap_masses, yield_mean_array - yield_sigma_array, yield_mean_array + yield_sigma_array)
            
            pdf.savefig()
            plt.close()

            for index in range(0, len(variable_map["sig_mean"])) : 
        
                plt.plot(ap_masses, variable_map["sig_pull"][index], 'o-', label=name_list[index])
                plt.xlabel("$A'$ mass hypothesis")
                plt.ylabel("Signal Yield Pull")
                plt.ylim([-2, 2])
                plt.legend()

            pdf.savefig()
            plt.close()
            '''


if __name__ == "__main__" : 
    main()

