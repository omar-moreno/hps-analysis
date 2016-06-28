
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import root_numpy as rnp

from matplotlib.backends.backend_pdf import PdfPages

class FeeAnalysis : 

    def __init__(self) : 

        # Use the 'Bayesian Methods for Hackers' style
        plt.style.use('bmh')

        # Set the default font size for plots and legends
        matplotlib.rcParams.update({'font.size': 10, 'legend.fontsize': 10})

        # Enable the use of LaTeX in titles
        plt.rc('text', usetex=True)

        self.cluster_energy = []
        self.cluster_size = []
        self.cluster_seed_energy = []
        self.cluster_seed_index_x = []
        self.cluster_seed_index_y = []
        self.cluster_time = []
        self.cluster_x = []
        self.cluster_y = []
        self.cluster_z = []
        self.has_track_match = []

    def process(self, root_file) : 

        print "[ FeeAnalysis ]: Processing file " + str(root_file) 
        results = rnp.root2array(root_file)

        self.cluster_energy = np.append(self.cluster_energy, results["cluster_energy"])
        print self.cluster_energy
        self.cluster_size = np.append(self.cluster_size, results["cluster_size"])
        self.cluster_seed_energy = np.append(self.cluster_seed_energy, results["cluster_seed_energy"])
        self.cluster_seed_index_x = np.append(self.cluster_seed_index_x, results["cluster_seed_index_x"])
        self.cluster_seed_index_y = np.append(self.cluster_seed_index_y, results["cluster_seed_index_y"])
        self.cluster_time = np.append(self.cluster_time, results["cluster_time"])
        self.cluster_x = np.append(self.cluster_x, results["cluster_x"])
        self.cluster_y = np.append(self.cluster_y, results["cluster_y"])
        self.cluster_z = np.append(self.cluster_z, results["cluster_z"])
        self.has_track_match = np.append(self.has_track_match, results["has_track_match"])

    def make_plots(self) : 
        with PdfPages("fee_analysis.pdf") as pdf :
            
            bins = np.linspace(0, 2.5, 250)
            plt.hist(self.cluster_energy, bins, alpha=0.8, histtype="stepfilled")
            plt.xlabel("Cluster energy (GeV)")
            pdf.savefig()
            plt.close()

            bins = np.linspace(0, 180, 360)
            plt.hist(self.cluster_time, bins, alpha=0.8, histtype="stepfilled")
            plt.xlabel("Cluster time (ns)")
            pdf.savefig()
            plt.close()
           
            bins_x = np.linspace(0, 2.5, 250)
            bins_y = np.linspace(0, 180, 360)
            plt.hist2d(self.cluster_energy, self.cluster_time, bins=(bins_x, bins_y))
            plt.xlabel("Cluster energy (GeV)")
            plt.ylabel("Cluster time (ns)")
            pdf.savefig()
            plt.close()
            
            bins_x = np.linspace(0, 2.5, 250)
            bins_y = np.linspace(0, 10, 10)
            plt.hist2d(self.cluster_energy, self.cluster_size, bins=(bins_x, bins_y))
            plt.xlabel("Cluster energy (GeV)")
            plt.ylabel("Cluster size")
            pdf.savefig()
            plt.close()
            
            bins_x = np.linspace(0, 180, 360)
            bins_y = np.linspace(0, 10, 10)
            plt.hist2d(self.cluster_time, self.cluster_size, bins=(bins_x, bins_y))
            plt.xlabel("Cluster time (ns)")
            plt.ylabel("Cluster size")
            pdf.savefig()
            plt.close()

            plt.hist2d(self.cluster_x, self.cluster_y)
            pdf.savefig()
            plt.close()
            
            bins_x = np.linspace(0, 2.5, 250)
            bins_y = np.linspace(0, 2.5, 250)
            plt.hist2d(self.cluster_energy, self.cluster_seed_energy, bins=(bins_x, bins_y))
            plt.xlabel("Cluster energy (GeV)")
            plt.ylabel("Cluster seed energy (GeV)")
            pdf.savefig()
            plt.close()
