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
    parser.add_argument("-l", "--input_list", help="ROOT file")
    args = parser.parse_args()

    if args.input_list is None : 
        print 'A list of ROOT files needs to be specified'
        sys.exit(2)

    file_list = open(args.input_list, 'r')
    results_array_list = []

    yield_mean_array_list = []
    yield_sigma_array_list = []
    yield_pull_array_list = []

    #labels = ['pol2', 'pol3', 'pol4', 'pol5', 'pol6']
    labels = ['pol3', 'pol5']

    for line in file_list : 
        
        print "Processing file: " + str(line.strip())

        results_rec = rnp.root2array(line.strip())
        results_array_list.append(rnp.rec2array(results_rec))
    
        yield_mean_array_list.append([])
        yield_sigma_array_list.append([])
        yield_pull_array_list.append([])

    ap_masses = np.unique(results_array_list[0][:,0])
    print "A' masses: " + str(ap_masses)

    with PdfPages("results.pdf") as pdf :

        for ap_mass in ap_masses :

            print "Processing A' mass: " + str(ap_mass)

            #fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, figsize=(15, 5))
            fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(10, 5))

            index = 0
            for results_array in results_array_list : 
            
                yield_array = results_array[:,1][results_array[:,0] == ap_mass]
                yield_error_array = results_array[:,2][results_array[:,0] == ap_mass]
            
                ax0.hist(yield_array, bins=50, alpha=0.8, histtype="stepfilled", normed=True, label=labels[index])
                ax0.set_xlabel('Background yield', fontsize=8)
                #ax0.set_xlabel('Signal yield', fontsize=8)
                ax0.set_ylabel('AU', fontsize=8)
                ax0.set_title("$A'$ mass hypothesis = " + str(ap_mass), fontsize=8)
                ax0.legend()

                mu, std = norm.fit(yield_array)
                yield_mean_array_list[index].append(mu)
                yield_sigma_array_list[index].append(std)
                xmin, xmax = ax0.get_xlim()
                x = np.linspace(xmin, xmax, 100)
                p = norm.pdf(x, mu, std)
                ax0.plot(x, p, linewidth=2)

                ax1.hist(yield_error_array, bins=50, alpha=0.8, histtype="stepfilled", normed=True, label=labels[index])
                ax1.set_xlabel('Background yield fit error', fontsize=8)
                #ax1.set_xlabel('Signal yield fit error', fontsize=8)
                ax1.set_ylabel('AU', fontsize=8)
                ax1.set_title("$A'$ mass hypothesis = " + str(ap_mass), fontsize=8)
                ax1.legend()

                #ax2.hist(yield_array/yield_error_array, bins=50, alpha=0.8, histtype="stepfilled", normed=True, label=labels[index])
                #ax2.set_xlabel('Pull', fontsize=8)
                #ax2.set_ylabel('AU', fontsize=8)
                #ax2.set_title("$A'$ mass hypothesis = " + str(ap_mass),fontsize=8)


                #mu, std = norm.fit(yield_array/yield_error_array)
                #yield_pull_array_list[index].append(mu)
                #xmin, xmax = ax2.get_xlim()
                #x = np.linspace(xmin, xmax, 100)
                #p = norm.pdf(x, mu, std)
                #ax2.plot(x, p, linewidth=2)

                index += 1

            pdf.savefig()
            plt.close()

        fig, (ax0, ax1) = plt.subplots(ncols=2)
            
        for index in range(0, len(yield_mean_array_list)) : 

            yield_mean_array = np.array(yield_mean_array_list[index])
            yield_sigma_array = np.array(yield_sigma_array_list[index])
            #yield_pull_array = np.array(yield_pull_array_list[index])

            ax0.plot(ap_masses, yield_mean_array, 'o-', label=labels[index])
            ax0.set_xlabel("$A'$ mass hypothesis")
            ax0.set_ylabel("Background yield")
            #ax0.set_ylabel("Signal yield")
            #ax0.fill_between(ap_masses, yield_mean_array - yield_sigma_array, yield_mean_array + yield_sigma_array)
            #ax0.set_ylim([-1100, 1100])
            ax0.legend()


            ax1.plot(ap_masses, yield_mean_array/yield_sigma_array, 'o-', label=labels[index])
            ax1.set_xlabel("$A'$ mass hypothesis")
            ax1.set_ylabel("Background yield pull")
            #ax1.set_ylabel("Signal yield pull")
            #ax1.set_ylim([-10, 10])
            ax1.legend()
            

        pdf.savefig()
        plt.close()

if __name__ == "__main__" : 
    main()

