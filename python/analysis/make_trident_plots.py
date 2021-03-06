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
    
    with PdfPages("trident_plots.pdf") as pdf : 
   
        #
        #  Ecal cluster positions
        #

        cluster_x_high_arr = results_rec["cluster_x_high"]
        cluster_x_low_arr = results_rec["cluster_x_low"]
        plt.hist2d(cluster_x_high_arr, cluster_x_low_arr, bins=300)
        
        pdf.savefig()
        plt.close()

        cluster_y_high_arr = results_rec["cluster_y_high"]
        cluster_y_low_arr = results_rec["cluster_y_low"]
        plt.hist2d(cluster_y_high_arr, cluster_y_low_arr, bins=300)
        
        pdf.savefig()
        plt.close()


        #
        # Chi2 Distributions
        #

        fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(ncols=2, nrows=2, figsize=(15, 15))
        
        chi2_bins = np.linspace(0, 60, 120)

        electron_chi2_arr = results_rec["electron_chi2"]
        electron_chi2_arr_cuts = electron_chi2_arr[electron_chi2_arr < 60]
        positron_chi2_arr = results_rec["positron_chi2"]
        positron_chi2_arr_cuts = positron_chi2_arr[positron_chi2_arr < 60]
       
        ax0.set_xlabel("$\chi^{2}$")
        ax1.set_xlabel("$\chi^{2}$")

        ax0.hist(electron_chi2_arr_cuts, chi2_bins, alpha=0.8, histtype="stepfilled",
                label="all")
        ax1.hist(positron_chi2_arr_cuts, chi2_bins, alpha=0.8, histtype="stepfilled",
                label="all")
        
        electron_chi2_arr_cuts = electron_chi2_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15)]
        positron_chi2_arr_cuts = positron_chi2_arr[
            (positron_chi2_arr < 15) & (electron_chi2_arr < 15)]

        ax0.hist(electron_chi2_arr_cuts, chi2_bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$")
        ax1.hist(positron_chi2_arr_cuts, chi2_bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$")
    
        ax0.legend(prop={'size':12})
        ax1.legend(prop={'size':12})

        ax2.hist2d(electron_chi2_arr, positron_chi2_arr, chi2_bins)
        
        ax3.hist2d(electron_chi2_arr_cuts, positron_chi2_arr_cuts, bins=100)

        pdf.savefig()
        plt.close()

        #
        # Track momentum distributions
        #
        
        fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(ncols=2, nrows=2, figsize=(15, 15))
        p_bins = np.linspace(0, 1.5, 200)
        
        e_p_arr = results_rec["electron_p"]
        p_p_arr = results_rec["positron_p"]
        e_p_arr_cuts = e_p_arr[(e_p_arr < 1.5) & (p_p_arr < 1.5)]
        p_p_arr_cuts = p_p_arr[(e_p_arr < 1.5) & (p_p_arr < 1.5)]
       
        ax0.set_xlabel("Electron p (GeV)")
        ax1.set_xlabel("Positron p (GeV)")

        ax0.hist(e_p_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled",
                label="all")
        ax1.hist(p_p_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled",
                label="all")
        ax2.hist2d(e_p_arr_cuts, p_p_arr_cuts, bins=300)

        e_p_arr_cuts = e_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 1.5) & (p_p_arr < 1.5)]
        p_p_arr_cuts = p_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 1.5) & (p_p_arr < 1.5)]
        ax0.hist(e_p_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$")
        ax1.hist(p_p_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$")

        e_p_arr_cuts = e_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85)]
        p_p_arr_cuts = p_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85)]
        ax0.hist(e_p_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$, $p_e^+ < 0.85$, $p_e^- < 0.85$")
        ax1.hist(p_p_arr_cuts, p_bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$, $p_e^+ < 0.85$, $p_e^- < 0.85$")

        ax3.hist2d(e_p_arr_cuts, p_p_arr_cuts, bins=300)
        ax3.set_xlabel("Electron p (GeV)")
        ax3.set_ylabel("Positron p (GeV)")

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

        vx_arr = results_rec["vx"]
        vx_arr_cuts = vx_arr[abs(vx_arr) < 0.5]
        ax0.hist(vx_arr_cuts, vx_bins, alpha=0.8, histtype="stepfilled",
                label="all")
        ax0.set_xlim(-0.5, 0.5)
        
        vy_bins = np.linspace(-0.1, 0.1, 200)
        vy_arr = results_rec["vy"]
        vy_arr_cuts = vy_arr[abs(vy_arr) < 0.1]
        ax1.hist(vy_arr_cuts, vy_bins, alpha=0.8, histtype="stepfilled",
                label="all")
        ax1.set_xlim(-0.1, 0.1)

        vx_arr_cuts = vx_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15)]
        ax0.hist(vx_arr_cuts, vx_bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$")

        vx_arr_cuts = vx_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85)]
        ax0.hist(vx_arr_cuts, vx_bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$, $p_e^+ < 0.85$, $p_e^- < 0.85$")

        vx_arr_cuts = vx_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)]
        ax0.hist(vx_arr_cuts, vx_bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$, $p_e^+ < 0.85$, $p_e^- < 0.85$, vertex cut")

        vy_arr_cuts = vy_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15)]
        ax1.hist(vy_arr_cuts, vy_bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$")

        vy_arr_cuts = vy_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85)]
        ax1.hist(vy_arr_cuts, vy_bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$, $p_e^+ < 0.85$, $p_e^- < 0.85$")

        vy_arr_cuts = vy_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)]
        ax1.hist(vy_arr_cuts, vy_bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$, $p_e^+ < 0.85$, $p_e^- < 0.85$, vertex cut")

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
        v0_p_arr = results_rec["v0_p"]
        plt.hist(v0_p_arr, bins, alpha=0.8, histtype="stepfilled",
                label="all")
        plt.xlabel("$V_{0} p$ (GeV)")
        
        v0_p_arr_cuts = v0_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$")

        v0_p_arr_cuts = v0_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$, $p_e^+ < 0.85$, $p_e^- < 0.85$")

        v0_p_arr_cuts = v0_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$, $p_e^+ < 0.85$, $p_e^- < 0.85$, vertex cut")

        v0_p_arr_cuts = v0_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1) &
            (v0_p_arr > 0.85)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$, $p_e^+ < 0.85$, $p_e^- < 0.85$, vertex cut, $p_{sum} > 0.85$")

        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=1, mode="expand", borderaxespad=0.)

        pdf.savefig()
        plt.close()

        bins = np.linspace(0, 0.1, 800)

        mass_arr = results_rec["invariant_mass"]
        mass_arr_cuts = mass_arr[mass_arr < 0.1]
        plt.hist(mass_arr, bins, alpha=0.8, histtype="stepfilled", 
                label="all")
        plt.xlabel("Invariant Mass M($e^-e^+$) (GeV)")
        plt.ylabel("Events/.2")
        plt.xlim(0, 0.1)
      
        mass_arr_cuts = mass_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15)]
        plt.hist(mass_arr_cuts, bins, alpha=0.8, histtype="stepfilled", 
                label="$\chi^2 < 15$")

        mass_arr_cuts = mass_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85)]
        plt.hist(mass_arr_cuts, bins, alpha=0.8, histtype="stepfilled", 
                label="$\chi^2 < 15$, $p_e^+ < 0.85$, $p_e^- < 0.85$")
        
        mass_arr_cuts = mass_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)]
        plt.hist(mass_arr_cuts, bins, alpha=0.8, histtype="stepfilled", 
                label="$\chi^2 < 15$, $p_e^+ < 0.85$, $p_e^- < 0.85$, vertex cut")

        mass_arr_cuts = mass_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1) &
            (v0_p_arr > 0.85)]
        plt.hist(mass_arr_cuts, bins, alpha=0.8, histtype="stepfilled",
                label="$\chi^2 < 15$, $p_e^+ < 0.85$, $p_e^- < 0.85$, vertex cut, $p_{sum} > 0.85$")

        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=1, mode="expand", borderaxespad=0.)
        
        pdf.savefig()
        plt.close()

        bins = np.linspace(0, 1.5, 500)
        v0_p_arr_cuts = v0_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled")

        v0_p_arr_cuts = v0_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85) &
            (abs(cluster_y_high_arr) > 30) & (abs(cluster_y_low_arr) > 30) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled")


        v0_p_arr_cuts = v0_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85) &
            (abs(cluster_y_high_arr) > 50) & (abs(cluster_y_low_arr) > 50) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled")
        
        v0_p_arr_cuts = v0_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85) &
            (abs(cluster_y_high_arr) > 70) & (abs(cluster_y_low_arr) > 70) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled")

        pdf.savefig()
        plt.close()


        labels = [
            "All",
            "Ecal cluster x position > 50 mm",
            "Ecal cluster x position > 100 mm",
            "Ecal cluster x position > 200 mm",
        ]

        bins = np.linspace(0, 1.5, 500)
        v0_p_arr_cuts = v0_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled", label=labels[0])
        plt.xlabel("V0_p (GeV)")

        v0_p_arr_cuts = v0_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85) &
            (abs(cluster_x_high_arr) > 50) & (abs(cluster_x_low_arr) > 50) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled", label=labels[1])

        v0_p_arr_cuts = v0_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85) &
            (abs(cluster_x_high_arr) > 100) & (abs(cluster_x_low_arr) > 100) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled", label=labels[2])

        v0_p_arr_cuts = v0_p_arr[
            (electron_chi2_arr < 15) & (positron_chi2_arr < 15) &
            (e_p_arr < 0.85) & (p_p_arr < 0.85) &
            (abs(cluster_x_high_arr) > 200) & (abs(cluster_x_low_arr) > 200) &
            (((vx_arr*vx_arr/(0.04)) + (vy_arr*vy_arr/(0.0025))) < 1)]
        plt.hist(v0_p_arr_cuts, bins, alpha=0.8, histtype="stepfilled", label=labels[3])
    
        plt.legend()
        

        pdf.savefig()
        plt.close()


    file = r.TFile("invariant_mass_final.root", "recreate")
    
    mass_histo = r.TH1F("invariant_mass", "invariant_mass", 800, 0., 0.1)
    for value in np.nditer(mass_arr_cuts) : 
        mass_histo.Fill(value)
        
    mass_histo.Write()

    f.Close()



if __name__ == "__main__" :
    main()
