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

            #rejection_efficiency = []
            #for index in xrange(0, 60) : 
            #    rejection_efficiency.append((len(


            bins = np.linspace(0, 180, 361)
            plt.hist(self.positron_cluster_time, bins, histtype="step")
            plt.hist(positron_cluster_time_acc, bins, histtype="step")
            plt.hist(positron_cluster_time_real, bins, histtype="step")
            pdf.savefig()
            plt.close()

            bins = np.linspace(0, 100, 101)
            plt.hist(self.electron_chi2, bins, histtype="step")
            plt.hist(electron_chi2_acc, bins, histtype="step")
            plt.hist(electron_chi2_real, bins, histtype="step")
            pdf.savefig()
            plt.close()
            
            bins = np.linspace(0, 1.5, 151)
            plt.hist(self.electron_p, bins, histtype="step")
            plt.hist(electron_p_acc, bins, histtype="step")
            plt.hist(electron_p_real, bins, histtype="step")
            pdf.savefig()
            plt.close()
            
            bins = np.linspace(0, 50, 51)
            plt.hist(self.v_chi2, bins, histtype="step")
            plt.hist(v_chi2_acc, bins, histtype="step")
            plt.hist(v_chi2_real, bins, histtype="step")
            pdf.savefig()
            plt.close()
