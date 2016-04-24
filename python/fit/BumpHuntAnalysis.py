import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import root_numpy as rnp

from  LogFormatterTex import LogFormatterTex
from matplotlib.ticker import LogFormatter
from matplotlib.backends.backend_pdf import PdfPages

class BumpHuntAnalysis : 

    def __init__(self) : 

        # Use the 'Bayesian Methods for Hackers' style
        plt.style.use('bmh')

        # Set the default font size for plots and legends
        matplotlib.rcParams.update({'font.size': 10, 'legend.fontsize': 10})


        # Enable the use of LaTeX in titles
        plt.rc('text', usetex=True)

        self.p_values = []
        self.mass = []
        self.upper_limits = []
        self.bkg_total = []
        self.bkg_window_size = []

    def process(self, root_file) : 

        print "[ BumpHuntAnalysis ]: Processing file " + str(root_file) 
        results = rnp.root2array(root_file)
        self.p_values.append(results["p_value"][0])
        self.mass.append(results["ap_mass"][0])
        self.upper_limits.append(results["upper_limit"][0])
        self.bkg_total.append(results["bkg_total"][0])
        self.bkg_window_size.append(results["bkg_window_size"][0])

    def make_plots(self) : 

        bkg_mev = np.divide(np.array(self.bkg_total), np.array(self.bkg_window_size)*1000)
        sig_bkg = np.divide(np.array(self.upper_limits)*.8, bkg_mev)
        epsilon = sig_bkg/(np.array(self.mass)*1000)
        epsilon = np.divide(epsilon*2, 137*3*math.pi*0.083)

        with PdfPages("bump_hunt_final_results.pdf") as pdf : 

            
            plt.plot(self.mass, self.p_values, marker='o', linestyle='-')
            plt.xlabel("$A'$ mass hypothesis (GeV)")
            plt.ylabel("p-value")
            plt.yscale('log')
            plt.xlim(0.02, 0.076)
            plt.ylim(0.05, 1)
            axes = plt.subplot(111)
            axes.yaxis.set_major_formatter(LogFormatterTex(labelOnlyBase=False))
            axes.yaxis.set_minor_formatter(LogFormatterTex(labelOnlyBase=False))
            plt.grid(True)
       
            pdf.savefig()
            plt.close()

            plt.plot(self.mass, self.upper_limits, marker='o', linestyle='-')
            plt.xlabel("$A'$ mass hypothesis (GeV)")
            plt.ylabel("Signal Upper Limit")
            plt.xlim(0.02, 0.076)
            pdf.savefig()
            plt.close()

            plt.plot(self.mass, bkg_mev, marker='o', linestyle='-')
            plt.xlabel("$A'$ mass hypothesis (GeV)")
            plt.ylabel("Background per MeV")
            plt.xlim(0.02, 0.076)
            pdf.savefig()
            plt.close()

            plt.plot(self.mass, sig_bkg, marker='o', linestyle='-')
            plt.xlabel("$A'$ mass hypothesis (GeV)")
            plt.ylabel("Signal/(Background per MeV)")
            plt.xlim(0.02, 0.076)
            plt.yscale('log')
            pdf.savefig()
            plt.close()

            plt.plot(self.mass, epsilon, marker='o', linestyle='-')
            plt.xlabel("$A'$ mass hypothesis (GeV)")
            plt.ylabel("$\epsilon^2$")
            plt.ylim(0.0000001, 0.001)
            plt.xlim(0.02, 0.076)
            plt.yscale('log')
            pdf.savefig()
            plt.close()
