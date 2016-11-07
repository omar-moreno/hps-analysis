from __future__ import division

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import root_numpy as rnp

from matplotlib.backends.backend_pdf import PdfPages
import ROOT as r

class TridentAnalysis : 

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

    def plt_sel_eff(self, ax, signal, min, max, step, invert) :

        cut = min
        sel_eff = []
        cuts    = []

        while cut <= max : 
            if invert : sig_integrated = len(signal[signal > cut])
            else : sig_integrated = len(signal[signal < cut])
            sel_eff.append((sig_integrated/len(signal))*100)
            cuts.append(cut)
            cut += step
        
        major_ticks = np.arange(0, 101, 10)                                              
        minor_ticks = np.arange(0, 101, 5)  

        ax.set_yticks(major_ticks)                                                       
        ax.set_yticks(minor_ticks, minor=True)  

        ax.grid(which='minor', alpha=0.2)                                                
        ax.grid(which='major', alpha=0.5)        
        
        ax.errorbar(cuts, sel_eff, marker='None', linestyle='-')
        ax.set_ylabel("Selection Efficiency (\%)")

    def plt_sig_bkg(self, ax, signal, bkg, min, max, step, invert) :

        cut = min
        sig_sqrt_sig_bkg = []
        cuts    = []

        while cut <= max : 
            if invert : 
                sig_integrated = len(signal[signal > cut])
                bkg_integrated = len(bkg[bkg > cut])
            else : 
                sig_integrated = len(signal[signal < cut])
                bkg_integrated = len(bkg[bkg < cut])
            sig_sqrt_sig_bkg.append(sig_integrated/math.sqrt(sig_integrated + bkg_integrated))
            cuts.append(cut)
            cut += step
        
        ax.errorbar(cuts, sig_sqrt_sig_bkg, marker='None', linestyle='-')
        ax.set_ylabel("$S/\sqrt(S+B)$")

        max_index = np.argmax(sig_sqrt_sig_bkg)
        print "Signal optimized at: " + str(cuts[max_index])

    def load_file(self, root_file) : 

        results = rnp.root2array(root_file, "results")

        cuts = results['electron_cluster_time'] != -9999
        self.electron_cluster_time = results['electron_cluster_time'][cuts]
        self.electron_chi2 = results['electron_chi2'][cuts]
        self.electron_p = results['electron_p'][cuts]
        
        self.positron_cluster_time = results['positron_cluster_time'][cuts]

        self.v_chi2 = results['v_chi2'][cuts]

        #self.mass = np.append(self.mass, results["invariant_mass"])
        #self.v0_p = np.append(self.v0_p, results["v0_p"])
        #self.e_p = np.append(self.e_p, results["electron_p"])
        #self.p_p = np.append(self.p_p, results["positron_p"])
        #self.e_py = np.append(self.e_py, results["electron_py"])
        #self.p_py = np.append(self.p_py, results["positron_py"])
    
    def process(self) :

        cluster_time_diff = self.electron_cluster_time - self.positron_cluster_time
        acc_cut = np.abs(cluster_time_diff) > 3
        sig_cut = np.abs(cluster_time_diff) < 1
        
        electron_cluster_time_acc = self.electron_cluster_time[acc_cut]
        positron_cluster_time_acc = self.positron_cluster_time[acc_cut]
        cluster_time_diff_acc = electron_cluster_time_acc - positron_cluster_time_acc
        electron_chi2_acc = self.electron_chi2[acc_cut]
        electron_p_acc = self.electron_p[acc_cut]
        v_chi2_acc = self.v_chi2[acc_cut]
        
        electron_cluster_time_real = self.electron_cluster_time[sig_cut]
        positron_cluster_time_real = self.positron_cluster_time[sig_cut]
        cluster_time_diff_real = electron_cluster_time_real - positron_cluster_time_real
        electron_chi2_real = self.electron_chi2[sig_cut]
        electron_p_real = self.electron_p[sig_cut]
        v_chi2_real = self.v_chi2[sig_cut]

        with PdfPages("v0_analysis.pdf") as pdf :
           
            bins = np.linspace(-10, 10, 201)
            plt.hist(cluster_time_diff, bins, histtype="step")
            plt.hist(cluster_time_diff_acc, bins, histtype="step")
            plt.hist(cluster_time_diff_real, bins, histtype="step")
            pdf.savefig()
            plt.close()

            bins = np.linspace(0, 180, 361)
            plt.hist(self.electron_cluster_time, bins, histtype="step")
            plt.hist(electron_cluster_time_acc, bins, histtype="step")
            plt.hist(electron_cluster_time_real, bins, histtype="step")
            pdf.savefig()
            plt.close()
           
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   

            self.plt_sel_eff(ax, electron_cluster_time_real, 20, 60, 0.0625, False)
            self.plt_sel_eff(ax, electron_cluster_time_real, 20, 60, 0.0625, True)
            ax.set_xlabel("Electron cluster time (ns)")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            self.plt_sig_bkg(ax, electron_cluster_time_real, electron_cluster_time_acc, 20, 60, 0.0625, False)
            self.plt_sig_bkg(ax, electron_cluster_time_real, electron_cluster_time_acc, 20, 60, 0.0625, True)
            ax.set_xlabel("Electron cluster time (ns)")
            pdf.savefig()
            plt.close()

            bins = np.linspace(0, 180, 361)
            plt.hist(self.positron_cluster_time, bins, histtype="step")
            plt.hist(positron_cluster_time_acc, bins, histtype="step")
            plt.hist(positron_cluster_time_real, bins, histtype="step")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_sel_eff(ax, positron_cluster_time_real, 20, 60, 0.0625, False)
            self.plt_sel_eff(ax, positron_cluster_time_real, 20, 60, 0.0625, True)
            ax.set_xlabel("Positron cluster time (ns)")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            self.plt_sig_bkg(ax, positron_cluster_time_real, positron_cluster_time_acc, 20, 60, 0.0625, False)
            self.plt_sig_bkg(ax, positron_cluster_time_real, positron_cluster_time_acc, 20, 60, 0.0625, True)
            ax.set_xlabel("Positron cluster time (ns)")
            pdf.savefig()
            plt.close()

            bins = np.linspace(0, 100, 101)
            plt.hist(self.electron_chi2, bins, histtype="step")
            plt.hist(electron_chi2_acc, bins, histtype="step")
            plt.hist(electron_chi2_real, bins, histtype="step")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_sel_eff(ax, electron_chi2_real, 1, 100, .5, False)
            ax.set_xlabel("Electron $\chi^2$")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            self.plt_sig_bkg(ax, electron_chi2_real, electron_chi2_acc, 1, 100, .5, False)
            ax.set_xlabel("Electron $\chi^2$")
            pdf.savefig()
            plt.close()


            bins = np.linspace(0, 1.5, 151)
            plt.hist(self.electron_p, bins, histtype="step")
            plt.hist(electron_p_acc, bins, histtype="step")
            plt.hist(electron_p_real, bins, histtype="step")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_sel_eff(ax, electron_p_real, .1, 1.4, .05, False)
            ax.set_xlabel("$p(e^-)$")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            self.plt_sig_bkg(ax, electron_p_real, electron_p_acc, .1, 1.4, .05, False)
            ax.set_xlabel("$p(e^-)$")
            pdf.savefig()
            plt.close()

            bins = np.linspace(0, 50, 51)
            plt.hist(self.v_chi2, bins, histtype="step")
            plt.hist(v_chi2_acc, bins, histtype="step")
            plt.hist(v_chi2_real, bins, histtype="step")
            pdf.savefig()
            plt.close()
