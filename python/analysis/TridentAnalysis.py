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

    def plt_sel_eff(self, ax, sig, bkg, min, max, step, invert) :

        cut = min
        sel_eff = []
        cuts    = []

        while cut <= max : 
            if invert : sig_integrated = len(sig[sig > cut])
            else : sig_integrated = len(sig[sig < cut])
            sel_eff.append((sig_integrated/len(sig))*100)
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

        return sel_eff

    def plt_rej_eff(self, ax, sig, bkg, min, max, step, invert) :

        cut = min
        rej_eff = []
        cuts    = []

        while cut <= max : 
            if invert : bkg_integrated = len(bkg[bkg > cut])
            else : bkg_integrated = len(bkg[bkg < cut])
            rej_eff.append((bkg_integrated/len(bkg))*100)
            cuts.append(cut)
            cut += step
        
        major_ticks = np.arange(0, 101, 10)                                              
        minor_ticks = np.arange(0, 101, 5)  

        ax.set_yticks(major_ticks)                                                       
        ax.set_yticks(minor_ticks, minor=True)  

        ax.grid(which='minor', alpha=0.2)                                                
        ax.grid(which='major', alpha=0.5)        
        
        ax.errorbar(cuts, rej_eff, marker='None', linestyle='-')
        ax.set_ylabel("Rejection Efficiency (\%)")

        return rej_eff


    def plt_roc_curve(self, ax, sel_eff, rej_eff) : 
        
        major_ticks = np.arange(0, 101, 10)                                              
        minor_ticks = np.arange(0, 101, 5)  

        ax.set_yticks(major_ticks)                                                       
        ax.set_yticks(minor_ticks, minor=True)  

        ax.grid(which='minor', alpha=0.2)                                                
        ax.grid(which='major', alpha=0.5)        
        
        ax.errorbar(sel_eff, 100 - np.array(rej_eff), marker='None', linestyle='-')
        ax.set_xlabel("Selection Efficiency (\%)")
        ax.set_ylabel("Rejection Efficiency (\%)")
        

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
        sig_cut_1pt0 = np.abs(cluster_time_diff) < 1
        sig_cut_0pt5 = np.abs(cluster_time_diff) < 0.5
        
        electron_cluster_time_acc = self.electron_cluster_time[acc_cut]
        positron_cluster_time_acc = self.positron_cluster_time[acc_cut]
        cluster_time_diff_acc = electron_cluster_time_acc - positron_cluster_time_acc
        electron_chi2_acc = self.electron_chi2[acc_cut]
        electron_p_acc = self.electron_p[acc_cut]
        v_chi2_acc = self.v_chi2[acc_cut]
        
        electron_cluster_time_real_1pt0 = self.electron_cluster_time[sig_cut_1pt0]
        positron_cluster_time_real_1pt0 = self.positron_cluster_time[sig_cut_1pt0]
        cluster_time_diff_real_1pt0 = electron_cluster_time_real_1pt0 - positron_cluster_time_real_1pt0
        electron_chi2_real_1pt0 = self.electron_chi2[sig_cut_1pt0]
        electron_p_real_1pt0 = self.electron_p[sig_cut_1pt0]
        v_chi2_real_1pt0 = self.v_chi2[sig_cut_1pt0]

        electron_cluster_time_real_0pt5 = self.electron_cluster_time[sig_cut_0pt5]
        positron_cluster_time_real_0pt5 = self.positron_cluster_time[sig_cut_0pt5]
        cluster_time_diff_real_0pt5 = electron_cluster_time_real_0pt5 - positron_cluster_time_real_0pt5
        electron_chi2_real_0pt5 = self.electron_chi2[sig_cut_0pt5]
        electron_p_real_0pt5 = self.electron_p[sig_cut_0pt5]
        v_chi2_real_0pt5 = self.v_chi2[sig_cut_0pt5]

        with PdfPages("trident_plots.pdf") as pdf :

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   

            bins = np.linspace(-10, 10, 201)
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            ax.hist(cluster_time_diff, bins, histtype="step", label="All")
            ax.hist(cluster_time_diff_acc, bins, histtype="step", label="$\Delta t >$ 3 ns")
            ax.hist(cluster_time_diff_real_1pt0, bins, histtype="step", label="$\Delta t <$ 1.0 ns")
            ax.hist(cluster_time_diff_real_0pt5, bins, histtype="step", label="$\Delta t <$ 0.5 ns")
            ax.set_xlabel("Electron cluster time - Positron cluster time (ns)")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))  

            bins = np.linspace(0, 180, 361)
            ax.hist(self.electron_cluster_time, bins, histtype="step", label="All")
            ax.hist(electron_cluster_time_acc, bins, histtype="step", label="$\Delta t >$ 3 ns")
            ax.hist(electron_cluster_time_real_1pt0, bins, histtype="step", label="$\Delta t >$ 1.0 ns")
            ax.hist(electron_cluster_time_real_0pt5, bins, histtype="step", label="$\Delta t >$ 0.5 ns")
            ax.set_xlabel("Electron cluster time (ns)")
            ax.legend()
            pdf.savefig()
            plt.close()
           
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   

            sel_eff_1pt0 = self.plt_sel_eff(ax, electron_cluster_time_real_1pt0, electron_cluster_time_acc, 0, 180, 2, False)
            print sel_eff_1pt0
            sel_eff_1pt0_invert = self.plt_sel_eff(ax, electron_cluster_time_real_1pt0, electron_cluster_time_acc, 0, 180, 0.0625, True)
            sel_eff_0pt5 = self.plt_sel_eff(ax, electron_cluster_time_real_0pt5, electron_cluster_time_acc, 0, 180, 0.0625, False)
            sel_eff_0pt5_invert = self.plt_sel_eff(ax, electron_cluster_time_real_0pt5, electron_cluster_time_acc, 0, 180, 0.0625, True)
            ax.set_xlabel("Electron cluster time (ns)")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   

            rej_eff_1pt0 = self.plt_rej_eff(ax, electron_cluster_time_real_1pt0, electron_cluster_time_acc, 0, 180, 2, False)
            print rej_eff_1pt0
            rej_eff_1pt0_invert = self.plt_rej_eff(ax, electron_cluster_time_real_1pt0, electron_cluster_time_acc, 0, 180, 0.0625, True)
            rej_eff_0pt5 = self.plt_rej_eff(ax, electron_cluster_time_real_0pt5, electron_cluster_time_acc, 0, 180, 0.0625, False)
            rej_eff_0pt5_invert = self.plt_rej_eff(ax, electron_cluster_time_real_0pt5, electron_cluster_time_acc, 0, 180, 0.0625, True)
            ax.set_xlabel("Electron cluster time (ns)")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_roc_curve(ax, sel_eff_1pt0, rej_eff_1pt0)
            self.plt_roc_curve(ax, sel_eff_0pt5, rej_eff_0pt5)
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_roc_curve(ax, sel_eff_1pt0_invert, rej_eff_1pt0_invert)
            self.plt_roc_curve(ax, sel_eff_0pt5_invert, rej_eff_0pt5_invert)
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_sig_bkg(ax, electron_cluster_time_real_1pt0, electron_cluster_time_acc, 20, 60, 0.0625, False)
            self.plt_sig_bkg(ax, electron_cluster_time_real_1pt0, electron_cluster_time_acc, 20, 60, 0.0625, True)
            self.plt_sig_bkg(ax, electron_cluster_time_real_0pt5, electron_cluster_time_acc, 20, 60, 0.0625, False)
            self.plt_sig_bkg(ax, electron_cluster_time_real_0pt5, electron_cluster_time_acc, 20, 60, 0.0625, True)
            ax.set_xlabel("Electron cluster time (ns)")
            pdf.savefig()
            plt.close()


            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   

            bins = np.linspace(0, 180, 361)
            ax.hist(self.positron_cluster_time, bins, histtype="step", label="All")
            ax.hist(positron_cluster_time_acc, bins, histtype="step", label="$\Delta t >$ 3 ns")
            ax.hist(positron_cluster_time_real_1pt0, bins, histtype="step", label="$\Delta t >$ 1.0 ns")
            ax.hist(positron_cluster_time_real_0pt5, bins, histtype="step", label="$\Delta t >$ 0.5 ns")
            ax.set_xlabel("Positron cluster time (ns)")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
           
            sel_eff_1pt0 = self.plt_sel_eff(ax, positron_cluster_time_real_1pt0, positron_cluster_time_acc, 20, 60, 0.0625, False)
            sel_eff_1pt0_invert = self.plt_sel_eff(ax, positron_cluster_time_real_1pt0, positron_cluster_time_acc, 20, 60, 0.0625, True)
            sel_eff_0pt5 = self.plt_sel_eff(ax, positron_cluster_time_real_0pt5, positron_cluster_time_acc, 20, 60, 0.0625, False)
            sel_eff_0pt5_invert = self.plt_sel_eff(ax, positron_cluster_time_real_0pt5, positron_cluster_time_acc, 20, 60, 0.0625, True)
            ax.set_xlabel("Positron cluster time (ns)")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            rej_eff_1pt0 = self.plt_rej_eff(ax, positron_cluster_time_real_1pt0, positron_cluster_time_acc, 20, 60, 0.0625, False)
            rej_eff_1pt0_invert = self.plt_rej_eff(ax, positron_cluster_time_real_1pt0, positron_cluster_time_acc, 20, 60, 0.0625, True)
            rej_eff_0pt5 = self.plt_rej_eff(ax, positron_cluster_time_real_0pt5, positron_cluster_time_acc, 20, 60, 0.0625, False)
            rej_eff_0pt5_invert = self.plt_rej_eff(ax, positron_cluster_time_real_0pt5, positron_cluster_time_acc, 20, 60, 0.0625, True)
            ax.set_xlabel("Positron cluster time (ns)")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_roc_curve(ax, sel_eff_1pt0, rej_eff_1pt0)
            self.plt_roc_curve(ax, sel_eff_0pt5, rej_eff_0pt5)
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_roc_curve(ax, sel_eff_1pt0_invert, rej_eff_1pt0_invert)
            self.plt_roc_curve(ax, sel_eff_0pt5_invert, rej_eff_0pt5_invert)
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_sig_bkg(ax, positron_cluster_time_real_1pt0, positron_cluster_time_acc, 20, 60, 0.0625, False)
            self.plt_sig_bkg(ax, positron_cluster_time_real_1pt0, positron_cluster_time_acc, 20, 60, 0.0625, True)
            self.plt_sig_bkg(ax, positron_cluster_time_real_0pt5, positron_cluster_time_acc, 20, 60, 0.0625, False)
            self.plt_sig_bkg(ax, positron_cluster_time_real_0pt5, positron_cluster_time_acc, 20, 60, 0.0625, True)
            ax.set_xlabel("Positron cluster time (ns)")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            bins = np.linspace(0, 100, 101)
            plt.hist(self.electron_chi2, bins, histtype="step", label="All")
            ax.hist(electron_chi2_acc, bins, histtype="step", label="$\Delta t >$ 3 ns")
            ax.hist(electron_chi2_real_1pt0, bins, histtype="step", label="$\Delta t >$ 1.0 ns")
            ax.hist(electron_chi2_real_0pt5, bins, histtype="step", label="$\Delta t >$ 0.5 ns")
            ax.set_xlabel("Electron $\chi^2$")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            sel_eff_1pt0 = self.plt_sel_eff(ax, electron_chi2_real_1pt0, electron_chi2_acc, 1, 100, .5, False)
            sel_eff_0pt5 = self.plt_sel_eff(ax, electron_chi2_real_0pt5, electron_chi2_acc, 1, 100, .5, False)
            ax.set_xlabel("Electron $\chi^2$")
            ax.legend()
            pdf.savefig()
            plt.close()

            rej_eff_1pt0 = self.plt_rej_eff(ax, electron_chi2_real_1pt0, electron_chi2_acc, 1, 100, .5, False)
            rej_eff_0pt5 = self.plt_rej_eff(ax, electron_chi2_real_0pt5, electron_chi2_acc, 1, 100, .5, False)
            ax.set_xlabel("Electron $\chi^2$")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_roc_curve(ax, sel_eff_1pt0, rej_eff_1pt0)
            self.plt_roc_curve(ax, sel_eff_0pt5, rej_eff_0pt5)
            pdf.savefig()
            plt.close()


            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_sig_bkg(ax, electron_chi2_real_1pt0, electron_chi2_acc, 1, 100, .5, False)
            self.plt_sig_bkg(ax, electron_chi2_real_0pt5, electron_chi2_acc, 1, 100, .5, False)
            ax.set_xlabel("Electron $\chi^2$")
            ax.legend()
            pdf.savefig()
            plt.close()

            #
            # Electron track momentum
            #

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            bins = np.linspace(0, 1.5, 151)
            ax.hist(self.electron_p, bins, histtype="step", label="All")
            ax.hist(electron_p_acc, bins, histtype="step", label="$\Delta t >$ 3 ns")
            ax.hist(electron_p_real_1pt0, bins, histtype="step", label="$\Delta t >$ 1.0 ns")
            ax.hist(electron_p_real_0pt5, bins, histtype="step", label="$\Delta t >$ 0.5 ns")
            ax.set_xlabel("$p(e^-)$ (GeV)")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            sel_eff_1pt0 = self.plt_sel_eff(ax, electron_p_real_1pt0, electron_p_acc, .1, 1.4, .05, False)
            sel_eff_0pt5 = self.plt_sel_eff(ax, electron_p_real_0pt5, electron_p_acc, .1, 1.4, .05, False)
            ax.set_xlabel("$p(e^-)$ (GeV)")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            rej_eff_1pt0 = self.plt_rej_eff(ax, electron_p_real_1pt0, electron_p_acc, .1, 1.4, .05, False)
            rej_eff_0pt5 = self.plt_rej_eff(ax, electron_p_real_0pt5, electron_p_acc, .1, 1.4, .05, False)
            ax.set_xlabel("$p(e^-)$ (GeV)")
            ax.legend()
            pdf.savefig()
            plt.close()
            
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_roc_curve(ax, sel_eff_1pt0, rej_eff_1pt0)
            self.plt_roc_curve(ax, sel_eff_0pt5, rej_eff_0pt5)
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_sig_bkg(ax, electron_p_real_1pt0, electron_p_acc, .1, 1.4, .05, False)
            self.plt_sig_bkg(ax, electron_p_real_0pt5, electron_p_acc, .1, 1.4, .05, False)
            ax.set_xlabel("$p(e^-)$ (GeV)")
            ax.legend()
            pdf.savefig()
            plt.close()
    
            #
            # TC vertex chi^2
            #

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            bins = np.linspace(0, 50, 51)
            ax.hist(self.v_chi2, bins, histtype="step", label="All")
            ax.hist(v_chi2_acc, bins, histtype="step", label="$\Delta t >$ 3 ns")
            ax.hist(v_chi2_real_1pt0, bins, histtype="step", label="$\Delta t >$ 1.0 ns")
            ax.hist(v_chi2_real_0pt5, bins, histtype="step", label="$\Delta t >$ 0.5 ns")
            ax.set_xlabel("Target-constrained vertex $\chi^2$")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            sel_eff_1pt0 = self.plt_sel_eff(ax, v_chi2_real_1pt0, v_chi2_acc, 1, 50, .5, False)
            sel_eff_0pt5 = self.plt_sel_eff(ax, v_chi2_real_0pt5, v_chi2_acc, 1, 50, .5, False)
            ax.set_xlabel("Target-constrained vertex $\chi^2$")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            rej_eff_1pt0 = self.plt_rej_eff(ax, v_chi2_real_1pt0, v_chi2_acc, 1, 50, .5, False)
            rej_eff_0pt5 = self.plt_rej_eff(ax, v_chi2_real_0pt5, v_chi2_acc, 1, 50, .5, False)
            ax.set_xlabel("Target-constrained vertex $\chi^2$")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_roc_curve(ax, sel_eff_1pt0, rej_eff_1pt0)
            self.plt_roc_curve(ax, sel_eff_0pt5, rej_eff_0pt5)
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_sig_bkg(ax, v_chi2_real_1pt0, v_chi2_acc, 1, 50, .5, False)
            self.plt_sig_bkg(ax, v_chi2_real_0pt5, v_chi2_acc, 1, 50, .5, False)
            ax.set_xlabel("Target-constrained vertex $\chi^2$")
            ax.legend()
            pdf.savefig()
            plt.close()
