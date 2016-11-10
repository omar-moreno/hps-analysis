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

        self.v_chi2 = results['v_chi2'][cuts]

        #self.mass = np.append(self.mass, results["invariant_mass"])
        #self.v0_p = np.append(self.v0_p, results["v0_p"])
        #self.e_p = np.append(self.e_p, results["electron_p"])
        #self.p_p = np.append(self.p_p, results["positron_p"])
        #self.e_py = np.append(self.e_py, results["electron_py"])
        #self.p_py = np.append(self.p_py, results["positron_py"])
    
    def process(self) :

        selection = self.electron_p < 0.8
        
        with PdfPages("trident_selection.pdf") as pdf :
        
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   
            
            bins = np.linspace(0, 1.5, 151)
            ax.hist(self.electron_p, bins, histtype="step", label="All")
            ax.hist(self.electron_p[selection], bins, histtype="step", label="$e(p) < 0.8")
