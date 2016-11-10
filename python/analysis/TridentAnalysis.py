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
        
        ax.errorbar(100 - np.array(rej_eff), sel_eff, marker='None', linestyle='-')
        ax.set_xlabel("Selection Efficiency (\%)")
        ax.set_ylabel("Rejection Efficiency (\%)")
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        

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
        self.electron_chi2 = results['electron_chi2'][cuts]
        self.electron_p = results['electron_p'][cuts]
        self.electron_time = results['electron_time'][cuts]
       
        self.electron_cluster_time = results['electron_cluster_time'][cuts]
        self.electron_cluster_x = results['electron_cluster_x'][cuts]
        self.electron_cluster_y = results['electron_cluster_y'][cuts]
        self.electron_cluster_z = results['electron_cluster_z'][cuts]

        self.positron_p = results['positron_p'][cuts]
        self.positron_time = results['positron_time'][cuts]
       
        self.positron_cluster_time = results['positron_cluster_time'][cuts]
        self.positron_cluster_x = results['positron_cluster_x'][cuts]
        self.positron_cluster_y = results['positron_cluster_y'][cuts]
        self.positron_cluster_z = results['positron_cluster_z'][cuts]

        self.top_chi2 = results['top_chi2'][cuts]
        self.top_p = results['top_p'][cuts]
        self.top_time = results['top_time'][cuts]
       
        self.top_cluster_time = results['top_cluster_time'][cuts]
        self.top_cluster_x = results['top_cluster_x'][cuts]
        self.top_cluster_y = results['top_cluster_y'][cuts]
        self.top_cluster_z = results['top_cluster_z'][cuts]

        self.bot_p = results['bot_p'][cuts]
        self.bot_time = results['bot_time'][cuts]
       
        self.bot_cluster_time = results['bot_cluster_time'][cuts]
        self.bot_cluster_x = results['bot_cluster_x'][cuts]
        self.bot_cluster_y = results['bot_cluster_y'][cuts]
        self.bot_cluster_z = results['bot_cluster_z'][cuts]


        self.v_chi2 = results['v_chi2'][cuts]

        #self.mass = np.append(self.mass, results["invariant_mass"])
        #self.v0_p = np.append(self.v0_p, results["v0_p"])
        #self.e_p = np.append(self.e_p, results["electron_p"])
        #self.p_p = np.append(self.p_p, results["positron_p"])
        #self.e_py = np.append(self.e_py, results["electron_py"])
        #self.p_py = np.append(self.p_py, results["positron_py"])
    
    def process(self) :
        
        # Define signal and background region
        cluster_time_diff = self.top_cluster_time - self.bot_cluster_time 
        acc_cut = np.abs(cluster_time_diff) > 3
        sig_cut = np.abs(cluster_time_diff) < 1
       
        top_cluster_time_acc = self.top_cluster_time[acc_cut]
        bot_cluster_time_acc = self.bot_cluster_time[acc_cut]
        cluster_time_diff_acc = top_cluster_time_acc - bot_cluster_time_acc
        
        electron_time_acc = self.electron_time[acc_cut]
        electron_track_cluster_time_diff_acc = electron_time_acc - top_cluster_time_acc
        positron_time_acc = self.positron_time[acc_cut]
        positron_track_cluster_time_diff_acc = positron_time_acc - bot_cluster_time_acc
        electron_chi2_acc = self.electron_chi2[acc_cut]
        electron_p_acc = self.electron_p[acc_cut]
        v_chi2_acc = self.v_chi2[acc_cut]
        
        top_cluster_time_sig = self.top_cluster_time[sig_cut]
        bot_cluster_time_sig = self.bot_cluster_time[sig_cut]
        cluster_time_diff_sig = top_cluster_time_sig - bot_cluster_time_sig

        electron_time_sig = self.electron_time[sig_cut]
        electron_track_cluster_time_diff_sig = electron_time_sig - top_cluster_time_sig
        
        positron_time_sig = self.positron_time[sig_cut]
        positron_track_cluster_time_diff_sig = positron_time_sig - bot_cluster_time_sig
        
        electron_chi2_sig = self.electron_chi2[sig_cut]
        electron_p_sig = self.electron_p[sig_cut]

        v_chi2_sig = self.v_chi2[sig_cut]

        with PdfPages("trident_plots.pdf") as pdf :

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   

            bins = np.linspace(-10, 10, 201)
            ax.hist(cluster_time_diff, bins, histtype="step", label="All")
            ax.hist(cluster_time_diff_acc, bins, histtype="step", label="$\Delta t >$ 3 ns")
            ax.hist(cluster_time_diff_sig, bins, histtype="step", label="$\Delta t <$ 1.0 ns")
            ax.set_xlabel("Top cluster time - Bottom cluster time (ns)")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))  

            bins_x = np.linspace(0, 1.5, 151)
            bins_y = np.linspace(0, 80, 161)
            ax.hist2d(self.top_p, self.top_cluster_time, bins=[bins_x, bins_y])
            pdf.savefig()
            plt.close()
           
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))  

            bins_x = np.linspace(0, 1.5, 151)
            bins_y = np.linspace(0, 80, 161)
            ax.hist2d(self.bot_p, self.bot_cluster_time, bins=[bins_x, bins_y])
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))  

            bins = np.linspace(0, 180, 361)
            ax.hist(self.top_cluster_time, bins, histtype="step", label="All")
            ax.hist(top_cluster_time_acc, bins, histtype="step", label="$\Delta t >$ 3 ns")
            ax.hist(top_cluster_time_sig, bins, histtype="step", label="$\Delta t >$ 1.0 ns")
            ax.set_xlabel("Top cluster time (ns)")
            ax.legend()
            pdf.savefig()
            plt.close()
           
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   

            sel_eff_1pt0 = self.plt_sel_eff(ax, top_cluster_time_sig, top_cluster_time_acc, 0, 180, 2, False)
            sel_eff_1pt0_invert = self.plt_sel_eff(ax, top_cluster_time_sig, top_cluster_time_acc, 0, 180, 0.0625, True)
            ax.set_xlabel("Top cluster time (ns)")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   

            rej_eff_1pt0 = self.plt_rej_eff(ax, top_cluster_time_sig, top_cluster_time_acc, 0, 180, 2, False)
            rej_eff_1pt0_invert = self.plt_rej_eff(ax, top_cluster_time_sig, top_cluster_time_acc, 0, 180, 0.0625, True)
            ax.set_xlabel("Top cluster time (ns)")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_roc_curve(ax, sel_eff_1pt0, rej_eff_1pt0)
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_roc_curve(ax, sel_eff_1pt0_invert, rej_eff_1pt0_invert)
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_sig_bkg(ax, top_cluster_time_sig, top_cluster_time_acc, 20, 60, 0.0625, False)
            self.plt_sig_bkg(ax, top_cluster_time_sig, top_cluster_time_acc, 20, 60, 0.0625, True)
            ax.set_xlabel("Electron cluster time (ns)")
            pdf.savefig()
            plt.close()



            #
            # Positron cluster time
            #

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   

            bins = np.linspace(0, 180, 361)
            ax.hist(self.bot_cluster_time, bins, histtype="step", label="All")
            ax.hist(bot_cluster_time_acc, bins, histtype="step", label="$\Delta t >$ 3 ns")
            ax.hist(bot_cluster_time_sig, bins, histtype="step", label="$\Delta t >$ 1.0 ns")
            ax.set_xlabel("Positron cluster time (ns)")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
           
            sel_eff_1pt0 = self.plt_sel_eff(ax, bot_cluster_time_sig, bot_cluster_time_acc, 20, 60, 0.0625, False)
            sel_eff_1pt0_invert = self.plt_sel_eff(ax, bot_cluster_time_sig, bot_cluster_time_acc, 20, 60, 0.0625, True)
            ax.set_xlabel("Positron cluster time (ns)")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            rej_eff_1pt0 = self.plt_rej_eff(ax, bot_cluster_time_sig, bot_cluster_time_acc, 20, 60, 0.0625, False)
            rej_eff_1pt0_invert = self.plt_rej_eff(ax, bot_cluster_time_sig, bot_cluster_time_acc, 20, 60, 0.0625, True)
            ax.set_xlabel("Positron cluster time (ns)")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_roc_curve(ax, sel_eff_1pt0, rej_eff_1pt0)
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_roc_curve(ax, sel_eff_1pt0_invert, rej_eff_1pt0_invert)
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_sig_bkg(ax, bot_cluster_time_sig, bot_cluster_time_acc, 20, 60, 0.0625, False)
            self.plt_sig_bkg(ax, bot_cluster_time_sig, bot_cluster_time_acc, 20, 60, 0.0625, True)
            ax.set_xlabel("Positron cluster time (ns)")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            bins = np.linspace(0, 100, 101)
            plt.hist(self.electron_chi2, bins, histtype="step", label="All")
            ax.hist(electron_chi2_acc, bins, histtype="step", label="$\Delta t >$ 3 ns")
            ax.hist(electron_chi2_sig, bins, histtype="step", label="$\Delta t >$ 1.0 ns")
            ax.set_xlabel("Electron $\chi^2$")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            sel_eff_1pt0 = self.plt_sel_eff(ax, electron_chi2_sig, electron_chi2_acc, 1, 100, .5, False)
            ax.set_xlabel("Electron $\chi^2$")
            ax.legend()
            pdf.savefig()
            plt.close()

            rej_eff_1pt0 = self.plt_rej_eff(ax, electron_chi2_sig, electron_chi2_acc, 1, 100, .5, False)
            ax.set_xlabel("Electron $\chi^2$")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_roc_curve(ax, sel_eff_1pt0, rej_eff_1pt0)
            pdf.savefig()
            plt.close()


            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_sig_bkg(ax, electron_chi2_sig, electron_chi2_acc, 1, 100, .5, False)
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
            ax.hist(electron_p_sig, bins, histtype="step", label="$\Delta t >$ 1.0 ns")
            ax.set_xlabel("$p(e^-)$ (GeV)")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            sel_eff_1pt0 = self.plt_sel_eff(ax, electron_p_sig, electron_p_acc, .1, 1.4, .05, False)
            ax.set_xlabel("$p(e^-)$ (GeV)")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            rej_eff_1pt0 = self.plt_rej_eff(ax, electron_p_sig, electron_p_acc, .1, 1.4, .05, False)
            ax.set_xlabel("$p(e^-)$ (GeV)")
            ax.legend()
            pdf.savefig()
            plt.close()
            
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_roc_curve(ax, sel_eff_1pt0, rej_eff_1pt0)
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_sig_bkg(ax, electron_p_sig, electron_p_acc, .1, 1.4, .05, False)
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
            ax.hist(v_chi2_sig, bins, histtype="step", label="$\Delta t >$ 1.0 ns")
            ax.set_xlabel("Target-constrained vertex $\chi^2$")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            sel_eff_1pt0 = self.plt_sel_eff(ax, v_chi2_sig, v_chi2_acc, 1, 50, .5, False)
            ax.set_xlabel("Target-constrained vertex $\chi^2$")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            rej_eff_1pt0 = self.plt_rej_eff(ax, v_chi2_sig, v_chi2_acc, 1, 50, .5, False)
            ax.set_xlabel("Target-constrained vertex $\chi^2$")
            ax.legend()
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_roc_curve(ax, sel_eff_1pt0, rej_eff_1pt0)
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            self.plt_sig_bkg(ax, v_chi2_sig, v_chi2_acc, 1, 50, .5, False)
            ax.set_xlabel("Target-constrained vertex $\chi^2$")
            ax.legend()
            pdf.savefig()
            plt.close()
