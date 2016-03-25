
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import root_numpy as rnp

from matplotlib.backends.backend_pdf import PdfPages

class BumpHuntToyAnalysis : 

    def __init__(self, root_file) : 

        # Use the 'Bayesian Methods for Hackers' style
        plt.style.use('bmh')

        # Set the default font size for plots and legends
        matplotlib.rcParams.update({'font.size': 16, 'legend.fontsize': 10})

        # Enable the use of LaTeX in titles
        plt.rc('text', usetex=True)

        self.results = rnp.root2array(root_file)

    def make_plots(self) : 
        
        # Create an array containing a unique list of histogram numbers
        histogram_numbers = np.unique(self.results["hist_n"])
        #print histogram_numbers

        hist_n = self.results["hist_n"]
        q0 = self.results["q0"]
        q0_max = []

        with PdfPages("bump_hunt_analysis.pdf") as pdf :

            # Find the largest q0 from all of the fits to a single histogram
            for histogram_number in histogram_numbers :

                q0_max.append(np.amax(q0[hist_n == histogram_number]))
            
            print q0_max
            bins = np.linspace(0, 10, 100)
            plt.hist(q0_max, bins, alpha=0.8, histtype="stepfilled")
            pdf.savefig()
            plt.close()




        
