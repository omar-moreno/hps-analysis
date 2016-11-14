from __future__ import division

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import root_numpy as rnp

from matplotlib.backends.backend_pdf import PdfPages
import ROOT as r

class TridentSelection : 

    def __init__(self) : 

        # Use the 'Bayesian Methods for Hackers' style
        plt.style.use('bmh')
        matplotlib.rcParams.update({'font.size': 12})
        matplotlib.rcParams['axes.facecolor'] = 'white'
        matplotlib.rcParams['legend.numpoints'] = 1        

        # Set the default font size for plots and legends
        matplotlib.rcParams.update({'font.size': 14, 'legend.fontsize': 14})

        # Enable the use of LaTeX in titles
        plt.rc('text', usetex=True)

    def load_file(self, root_file) : 

        results = rnp.root2array(root_file, "results")

        cuts = results['electron_p'] != -9999
        self.electron_chi2 = results['electron_chi2'][cuts]
        self.electron_p    = results['electron_p'][cuts]
        self.electron_time = results['electron_time'][cuts]
       
        self.positron_p    = results['positron_p'][cuts]
        self.positron_chi2 = results['positron_chi2'][cuts]
        self.positron_time = results['positron_time'][cuts]

        self.p_sum = np.array(self.electron_p) + np.array(self.positron_p)
       
        self.v_chi2 = results['v_chi2'][cuts]
        self.mass = results['invariant_mass'][cuts]
        self.v0_p =results["v0_p"][cuts]
    
    def process(self) :

        chi2_cut      = (self.electron_chi2 < 40) & (self.positron_chi2 < 40)
        fee_cut       = self.electron_p < 0.8
        p_sum_cut     = self.v0_p < 1.2
        radiative_cut = self.v0_p > 0.8
        selection     = chi2_cut & fee_cut & p_sum_cut & radiative_cut
        
        with PdfPages("trident_selection.pdf") as pdf :
       
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))  
            bins = np.linspace(0, 100, 101)
            ax.hist(self.electron_chi2, bins, histtype="step")
            ax.hist(self.electron_chi2[chi2_cut], bins, histtype="step")
            ax.set_xlabel("$e^-$ $\chi^2$")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))  
            bins = np.linspace(0, 100, 101)
            ax.hist(self.positron_chi2, bins, histtype="step")
            ax.hist(self.positron_chi2[chi2_cut], bins, histtype="step")
            ax.set_xlabel("$e^+$ $\chi^2$")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            bins = np.linspace(0, 1.5, 151)
            ax.hist(self.electron_p, bins, histtype="step", label="All")
            ax.hist(self.electron_p[fee_cut], bins, histtype="step", label="$e(p) < 0.8")
            pdf.savefig()
            plt.close()

            
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5)) 
            bins = np.linspace(0, 1.5, 151)
            ax.hist(self.p_sum, bins, histtype="step", label="All")
            ax.hist(self.p_sum[selection], bins, histtype="step", label="$e(p) < 0.8")
            pdf.savefig()
            plt.close()
            
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5)) 
            bins = np.linspace(0, .1, 2001)
            ax.hist(self.mass, bins, histtype="step")
            ax.hist(self.mass[selection], bins, histtype="step")
            pdf.savefig()
            plt.close()

            file = r.TFile("invariant_mass_final_selection.root", "recreate")
            mass_histo = r.TH1F("invariant_mass", "invariant_mass", 2000, 0., 0.1)
            for value in np.nditer(self.mass[selection]) : 
                mass_histo.Fill(value)
        
            mass_histo.Write()
            file.Close()

