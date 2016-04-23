#!/usr/bin/python

import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import root_numpy as rnp
import ROOT as r
import sys

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Ellipse
from rootpy.io import root_open
from scipy.stats import norm

def main() : 
   
    # Use the Bayesian Methods for Hackers design
    plt.style.use('bmh')
    matplotlib.rcParams.update({'font.size': 16, 'legend.fontsize': 10})
    plt.rc('text', usetex=True)

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-f", "--root_file", help="ROOT file to process.")
    args = parser.parse_args()

    # If a list of input files has not been specified, warn the user and exit 
    # the app.
    if args.root_file is None : 
        print 'Failed to specify a ROOT file.'
        sys.exit(2)

    results_rec = rnp.root2array(args.root_file)

#   features = rnp.list_branches(args.root_file)
#   print '\n'.join(str(feature) for feature in features)
    
    with PdfPages("moller_plots.pdf") as pdf : 

        n_tracks_arr = results_rec["n_tracks"]
        cluster_x_high_arr = results_rec["cluster_x_high"]
        cluster_x_low_arr = results_rec["cluster_x_low"]
        cluster_y_high_arr = results_rec["cluster_y_high"]
        cluster_y_low_arr = results_rec["cluster_y_low"]
        e_h_arr = results_rec["electron_high_p"]
        e_l_arr = results_rec["electron_low_p"]
        cluster_e_h = results_rec["cluster_energy_high"]
        cluster_e_l = results_rec["cluster_energy_low"]
        electron_high_chi2_arr = results_rec["electron_high_chi2"]
        electron_low_chi2_arr = results_rec["electron_low_chi2"]
        v_chi2_arr = results_rec["v_chi2"]
        vx_arr = results_rec["vx"]
        vy_arr = results_rec["vy"]
        v0_p_arr = results_rec["v0_p"]
        mass_arr = results_rec["invariant_mass"]
   
        #
        #  Ecal cluster x position
        #
        fig, (ax0, ax1) = plt.subplots(ncols=2, nrows=1, figsize=(20, 10))

        ax0.hist2d(cluster_x_high_arr, cluster_x_low_arr, bins=300)
        
        # Require both clusters to be on the electron side
        cluster_x_high_arr_cut = cluster_x_high_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0)
        ]
        cluster_x_low_arr_cut = cluster_x_low_arr[
            (cluster_x_low_arr < 0) & (cluster_x_high_arr < 0)
        ]

        ax1.hist2d(cluster_x_high_arr_cut, cluster_x_low_arr_cut, bins=300)

        pdf.savefig()
        plt.close()

        #
        # Ecal cluster x sum and diff
        #

        fig, (ax0, ax1) = plt.subplots(ncols=2, nrows=1, figsize=(20, 10))

        cluster_x_sum_arr_cut = cluster_x_high_arr_cut + cluster_x_low_arr_cut
        ax0.hist(cluster_x_sum_arr_cut, bins=300, alpha=0.8, histtype="stepfilled")

        cluster_x_diff_arr_cut = cluster_x_high_arr_cut - cluster_x_low_arr_cut
        ax1.hist(cluster_x_diff_arr_cut, bins=300, alpha=0.8, histtype="stepfilled")

        pdf.savefig()
        plt.close()

        #
        # Ecal cluster y
        #

        plt.hist2d(cluster_y_high_arr, cluster_y_low_arr, bins=300)
        
        pdf.savefig()
        plt.close()
    
        # 
        # E/p 
        #

        e_h_arr_cuts = e_h_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) 
        ]
        e_l_arr_cuts = e_l_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) 
        ]
        cluster_e_h_cuts = cluster_e_h[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) 
        ]
        cluster_e_l_cuts = cluster_e_l[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) 
        ]
        
        ep_h_arr_cuts = cluster_e_h_cuts/e_h_arr_cuts
        ep_l_arr_cuts = cluster_e_l_cuts/e_l_arr_cuts
        
        fig, (ax0, ax1) = plt.subplots(ncols=2, nrows=1, figsize=(20, 10))
        
        ep_bins = np.linspace(0, 2, 200)

        ax0.hist(ep_h_arr_cuts, ep_bins, alpha=0.8, histtype="stepfilled")
        ax0.set_xlim(0, 2)
        ax1.hist(ep_l_arr_cuts, ep_bins, alpha=0.8, histtype="stepfilled")
        ax1.set_xlim(0, 2)


        e_h_arr_cuts = e_h_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1) &
            (v0_p_arr > 1.056*.8) & (v0_p_arr < 1.2) & 
            (v_chi2_arr < 10)
        ]
        e_l_arr_cuts = e_l_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1) &
            (v0_p_arr > 1.056*.8) & (v0_p_arr < 1.2) & 
            (v_chi2_arr < 10)
        ]
        cluster_e_h_cuts = cluster_e_h[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1) &
            (v0_p_arr > 1.056*.8) & (v0_p_arr < 1.2) & 
            (v_chi2_arr < 10)
        ]
        cluster_e_l_cuts = cluster_e_l[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1) &
            (v0_p_arr > 1.056*.8) & (v0_p_arr < 1.2) & 
            (v_chi2_arr < 10)
        ]
        
        ep_h_arr_cuts = cluster_e_h_cuts/e_h_arr_cuts
        ep_l_arr_cuts = cluster_e_l_cuts/e_l_arr_cuts

        ax0.hist(ep_h_arr_cuts, ep_bins, alpha=0.8, histtype="stepfilled")
        ax0.set_xlim(0, 2)
        ax1.hist(ep_l_arr_cuts, ep_bins, alpha=0.8, histtype="stepfilled")
        ax1.set_xlim(0, 2)

        pdf.savefig()
        plt.close()

        #
        # Chi2 Distributions
        #
        fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(ncols=2, nrows=2, figsize=(15, 15))
        
        chi2_bins = np.linspace(0, 60, 120)

        electron_high_chi2_arr_cuts = electron_high_chi2_arr[electron_high_chi2_arr < 60]
        electron_low_chi2_arr_cuts = electron_low_chi2_arr[electron_low_chi2_arr < 60]
       
        ax0.set_xlabel("Electron Track $\chi^{2}$ - Highest Energy")
        ax1.set_xlabel("Electron Track $\chi^{2}$ - Lowest Energy")

        ax0.hist(electron_high_chi2_arr_cuts, chi2_bins, alpha=0.8, histtype="stepfilled",label="all")
        ax1.hist(electron_low_chi2_arr_cuts, chi2_bins, alpha=0.8, histtype="stepfilled", label="all")

        electron_high_chi2_arr_cuts = electron_high_chi2_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0)
        ]
        electron_low_chi2_arr_cuts = electron_low_chi2_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0)
        ]
        ax0.hist(electron_high_chi2_arr_cuts, chi2_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$")
        ax1.hist(electron_low_chi2_arr_cuts, chi2_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$")

        electron_high_chi2_arr_cuts = electron_high_chi2_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_high_chi2_arr < 15) & (electron_low_chi2_arr < 15)
        ]
        electron_low_chi2_arr_cuts = electron_low_chi2_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15)
        ]

        ax0.hist(electron_high_chi2_arr_cuts, chi2_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$")
        ax1.hist(electron_low_chi2_arr_cuts, chi2_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$")
    
        ax0.legend(prop={'size':12})
        ax1.legend(prop={'size':12})

        ax2.hist2d(electron_high_chi2_arr, electron_low_chi2_arr, chi2_bins)
        
        ax3.hist2d(electron_high_chi2_arr_cuts, electron_low_chi2_arr_cuts, bins=100)

        pdf.savefig()
        plt.close()

        #
        # Track momentum distributions
        #

        fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(ncols=2, nrows=2, figsize=(15, 15))
        p_bins = np.linspace(0, 1.5, 200)
        
        e_h_arr_cuts = e_h_arr[(e_h_arr < 1.5) & (e_l_arr < 1.5)]
        e_l_arr_cuts = e_l_arr[(e_h_arr < 1.5) & (e_l_arr < 1.5)]
       
        ax0.set_xlabel("Electron $p$ - Highest Energy (GeV)")
        ax1.set_xlabel("Electron $p$ - Lowest Energy (GeV)")

        ax0.hist(e_h_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled", label="all")
        ax1.hist(e_l_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled", label="all")
        ax2.hist2d(e_h_arr_cuts, e_l_arr_cuts, bins=300)

        e_h_arr_cuts = e_h_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (e_h_arr < 1.5) & (e_l_arr < 1.5)
        ]
        e_l_arr_cuts = e_l_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (e_h_arr < 1.5) & (e_l_arr < 1.5)
        ]
        ax0.hist(e_h_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$")
        ax1.hist(e_l_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$")

        e_h_arr_cuts = e_h_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15)
        ]
        e_l_arr_cuts = e_l_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15)
        ]
        ax0.hist(e_h_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$")
        ax1.hist(e_l_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$")

        e_h_arr_cuts = e_h_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70)
        ]
        e_l_arr_cuts = e_l_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70)
        ]
        ax0.hist(e_h_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$ \& $p(e) < 0.70$")
        ax1.hist(e_l_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$ \& $p(e) < 0.70$")

        ax3.hist2d(e_h_arr_cuts, e_l_arr_cuts, bins=300)
        ax3.set_xlabel("Electron $p$ - Highest Energy (GeV)")
        ax3.set_ylabel("Electron $p$ - Lowest Energy (GeV)")

        ax0.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=1, mode="expand", borderaxespad=0.)
        ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=1, mode="expand", borderaxespad=0.)

        pdf.savefig()
        plt.close()

        #
        # Vertex distributions
        #

        fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(ncols=2, nrows=2, figsize=(15, 15))

        ax0.set_xlabel("$V_{x}$ (mm)")
        ax1.set_xlabel("$V_{y}$ (mm)")

        vx_bins = np.linspace(-0.5, 0.5, 200)

        vx_arr_cuts = vx_arr[abs(vx_arr) < 0.5]
        ax0.hist(vx_arr_cuts, vx_bins, alpha=0.8, histtype="stepfilled", label="all")
        ax0.set_xlim(-0.5, 0.5)
        
        vy_bins = np.linspace(-0.1, 0.1, 200)
        vy_arr_cuts = vy_arr[abs(vy_arr) < 0.1]
        ax1.hist(vy_arr_cuts, vy_bins, alpha=0.8, histtype="stepfilled", label="all")
        ax1.set_xlim(-0.1, 0.1)

        vx_arr_cuts = vx_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0)
        ]
        ax0.hist(vx_arr_cuts, vx_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$")

        vx_arr_cuts = vx_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15)
        ]
        ax0.hist(vx_arr_cuts, vx_bins, alpha=0.8, histtype="stepfilled",
                 label="Ecal cluster x $< 0$ \& $\chi^2 < 15$")

        vx_arr_cuts = vx_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70)
        ]
        ax0.hist(vx_arr_cuts, vx_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$ \& $p(e) < 0.70$")

        vx_arr_cuts = vx_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)
        ]
        ax0.hist(vx_arr_cuts, vx_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$ \& $p(e) < 0.70$ \& vertex")

        vx_arr_cuts = vx_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1) &
            (v_chi2_arr < 10)
        ]
        ax0.hist(vx_arr_cuts, vx_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$ \& $p(e) < 0.70$ \& vertex")

        vy_arr_cuts = vy_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0)
        ]
        ax1.hist(vy_arr_cuts, vy_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$")

        vy_arr_cuts = vy_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15)
        ]
        ax1.hist(vy_arr_cuts, vy_bins, alpha=0.8, histtype="stepfilled",
                 label="Ecal cluster x $< 0$ \& $\chi^2 < 15$")

        vy_arr_cuts = vy_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70)
        ]
        ax1.hist(vy_arr_cuts, vy_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$ \& $p(e) < 0.70$")

        vy_arr_cuts = vy_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)
        ]
        ax1.hist(vy_arr_cuts, vy_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$ \& $p(e) < 0.70$ \& vertex")

        vy_arr_cuts = vy_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1) &
            (v_chi2_arr < 10)
        ]
        ax1.hist(vy_arr_cuts, vy_bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$ \& $p(e) < 0.70$ \& vertex")

        ax2.hist2d(vx_arr_cuts, vy_arr_cuts, bins=300)
 
        ax3.hist2d(results_rec["vx"], results_rec["vy"], bins=500)
        ax3.set_xlabel("$V_{x}$ (mm)")
        ax3.set_ylabel("$V_{y}$ (mm)")
        ax3.set_xlim(-0.5, 0.5)
        ax3.set_ylim(-0.2, 0.2)
    
        ellipse = Ellipse(xy=(0, 0), width=0.4, height=0.1, fc='None', edgecolor='r')
        ax3.add_artist(ellipse);
        
        ax0.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=1, mode="expand", borderaxespad=0.)
        ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=1, mode="expand", borderaxespad=0.)
        
        pdf.savefig()
        plt.close()

        #
        # V0 p
        #

        bins = np.linspace(0, 1.5, 500)
        plt.hist(v0_p_arr, bins, alpha=0.8, histtype="stepfilled", label="all")
        plt.xlabel("$V_{0} p$ (GeV)")

        v0_p_arr_cuts = v0_p_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0)
        ]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$")

        v0_p_arr_cuts = v0_p_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_high_chi2_arr < 15) & (electron_low_chi2_arr < 15)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                 label="Ecal cluster x $< 0$ \& $\chi^2 < 15$")

        v0_p_arr_cuts = v0_p_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70)
        ]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$ \& $p(e) < 0.70$")

        v0_p_arr_cuts = v0_p_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)
        ]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$ \& $p(e) < 0.70$ \& vertex")

        v0_p_arr_cuts = v0_p_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1) &
            (v0_p_arr > 1.056*.8) & (v0_p_arr < 1.2)
        ]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$ \& $p(e) < 0.70$ \& vertex \& $p_{sum}$")

        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=1, mode="expand", borderaxespad=0.)

        pdf.savefig()
        plt.close()

        e_h_arr_cuts = e_h_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (v0_p_arr > 1.056*.8) & (v0_p_arr < 1.2)
        ]
        e_l_arr_cuts = e_l_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (v0_p_arr > 1.056*.8) & (v0_p_arr < 1.2)
        ]
        plt.hist2d(e_h_arr_cuts, e_l_arr_cuts, bins=300)
        plt.xlabel("Electron $p$ - Highest Energy (GeV)")
        plt.ylabel("Electron $p$ - Lowest Energy (GeV)")

        pdf.savefig()
        plt.close()
    
        #
        # Invariant Mass
        #

        bins = np.linspace(0, 0.1, 800)

        mass_arr_cuts = mass_arr[mass_arr < 0.1]
        plt.hist(mass_arr, bins, alpha=0.8, histtype="stepfilled", label="all")
        plt.xlabel("Invariant Mass M($e^-e^+$) (GeV)")
        plt.ylabel("Events/.2")
        plt.xlim(0, 0.1)

        mass_arr_cuts = mass_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0)
        ]
        plt.hist(mass_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$")

        mass_arr_cuts = mass_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_high_chi2_arr < 15) & (electron_low_chi2_arr < 15)]
        plt.hist(mass_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                 label="Ecal cluster x $< 0$ \& $\chi^2 < 15$")

        mass_arr_cuts = mass_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70)
        ]
        plt.hist(mass_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$ \& $p(e) < 0.70$")

        mass_arr_cuts = mass_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)
        ]
        plt.hist(mass_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$ \& $p(e) < 0.70$ \& vertex")

        mass_arr_cuts = mass_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1) &
            (v0_p_arr > 1.056*.8) & (v0_p_arr < 1.2) & 
            (v_chi2_arr < 10)
        ]
        plt.hist(mass_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$ \& $p(e) < 0.70$ \& vertex \& $p_{sum}$")

        '''
        mass_arr_cuts = mass_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1) &
            (v0_p_arr > 1.056*.8) & (v0_p_arr < 1.2) &
            (n_tracks_arr > 2)
        ]
        plt.hist(mass_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                label="Ecal cluster x $< 0$ \& $\chi^2 < 15$ \& $p(e) < 0.70$ \& vertex \& $p_{sum}$")
        '''
        
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                     ncol=1, mode="expand", borderaxespad=0.)

        pdf.savefig()
        plt.close()

        cluster_x_high_arr_cut = cluster_x_high_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1) &
            (v0_p_arr > 1.056*.8) & (v0_p_arr < 1.2) &
            (v_chi2_arr < 10)
        ]
        cluster_x_low_arr_cut = cluster_x_low_arr[
            (cluster_x_high_arr < 0) & (cluster_x_low_arr < 0) &
            (electron_low_chi2_arr < 15) & (electron_high_chi2_arr < 15) &
            (e_h_arr < 0.70) & (e_l_arr < 0.70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1) &
            (v0_p_arr > 1.056*.8) & (v0_p_arr < 1.2) &
            (v_chi2_arr < 10)
        ]            

        plt.hist2d(cluster_x_high_arr_cut, cluster_x_low_arr_cut, bins=300)

        pdf.savefig()
        plt.close()
       
        fig, (ax0, ax1) = plt.subplots(ncols=2, nrows=1, figsize=(20, 10))

        cluster_x_sum_arr_cut = cluster_x_high_arr_cut + cluster_x_low_arr_cut
        ax0.hist(cluster_x_sum_arr_cut, bins=300, alpha=0.8, histtype="stepfilled")

        cluster_x_diff_arr_cut = cluster_x_high_arr_cut - cluster_x_low_arr_cut
        ax1.hist(cluster_x_diff_arr_cut, bins=300, alpha=0.8, histtype="stepfilled")

        pdf.savefig()
        plt.close()

        plt.hist2d(mass_arr_cuts, cluster_x_sum_arr_cut, bins=500)
        
        pdf.savefig()
        plt.close()

        plt.hist2d(mass_arr_cuts, cluster_x_diff_arr_cut, bins=500)
        
        pdf.savefig()
        plt.close()
 
        fig, (ax0, ax1) = plt.subplots(ncols=2, nrows=1, figsize=(20, 10))
       
        mass_arr_cuts = mass_arr_cuts[
            (ep_h_arr_cuts < 2) & (ep_l_arr_cuts < 2) & (ep_l_arr_cuts > 0)
        ]
        ep_h_arr_cutss = ep_h_arr_cuts[
            (ep_h_arr_cuts < 2) & (ep_l_arr_cuts < 2) & (ep_l_arr_cuts > 0)
        ]
        ep_l_arr_cutss = ep_l_arr_cuts[
            (ep_h_arr_cuts < 2) & (ep_l_arr_cuts < 2) & (ep_l_arr_cuts > 0)
        ]
        ax0.hist2d(mass_arr_cuts, ep_h_arr_cutss, bins=500)
        ax0.set_ylim(0,2)
        ax1.hist2d(mass_arr_cuts, ep_l_arr_cutss, bins=500)
        ax1.set_ylim(0,2)

        pdf.savefig()
        plt.close()

    root_file = r.TFile("moller_invariant_mass.root", "recreate")
    
    mass_histo = r.TH1F("invariant_mass", "invariant_mass", 800, 0., 0.1)
    for value in np.nditer(mass_arr_cuts) : 
        mass_histo.Fill(value)
        
    mass_histo.Write()

    root_file.Close()

if __name__ == "__main__" :
    main()
