
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
        print histogram_numbers

        #with PdfPages("bump_hunt_analysis.pdf") as pdf :

            # Find the largest




        
