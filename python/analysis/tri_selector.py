#!/usr/bin/python

from __future__ import division

import argparse
import sys
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import root_numpy as rnp
import ROOT as r

from matplotlib.backends.backend_pdf import PdfPages

def main() :

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-l", "--file_list", help="List of ROOT files to process.")
    args = parser.parse_args()

    if not args.file_list:
        print 'A list of ROOT files to process needs to be specified.'
        sys.exit(2)

    # Open the file containing the list of files to process
    root_file_list = None
    try:
        root_file_list = open(args.file_list, 'r')
    except IOError: 
        print "Unable to open file %s" % args.file_list
        sys.exit(2)

    root_files = []
    for line in root_file_list: 
        root_files.append(line.strip())

    rec = rnp.root2array(root_files, 'results')

    apply_tri_selection(rec)

def apply_tri_selection(rec):

    # Use the 'Bayesian Methods for Hackers' style
    plt.style.use('bmh')
    matplotlib.rcParams.update({'font.size': 12})
    matplotlib.rcParams['axes.facecolor'] = 'white'
    matplotlib.rcParams['legend.numpoints'] = 1        

    # Set the default font size for plots and legends    
    matplotlib.rcParams.update({'font.size': 16, 'legend.fontsize': 9})

    # Enable the use of LaTeX in titles
    plt.rc('text', usetex=True)

    electron_p      = rec['electron_p']
    electron_px     = rec['electron_px']
    electron_py     = rec['electron_py']
    electron_chi2   = rec['electron_chi2']
    electron_has_l1 = rec['electron_has_l1']
    electron_pt = np.sqrt(np.power(electron_px, 2) + np.power(electron_py, 2))
    
    positron_p      = rec['positron_p']
    positron_px     = rec['positron_px']
    positron_py     = rec['positron_py']
    positron_chi2   = rec['positron_chi2']
    positron_d0     = rec['positron_d0']
    positron_has_l1 = rec['positron_has_l1']
    positron_has_l2 = rec['positron_has_l2']

    positron_pt = np.sqrt(np.power(positron_px, 2) + np.power(positron_py, 2))

    top_cluster_time  = rec['top_cluster_time']
    bot_cluster_time  = rec['bot_cluster_time']
    cluster_time_diff = top_cluster_time - bot_cluster_time
    
    top_time = rec['top_time']
    
    bot_time = rec['bot_time']

    mass = rec['invariant_mass']
    v0_p = rec["v0_p"]


    #
    # Define cuts
    #
    # Accidentals
    radiative_cut = v0_p > 0.8*1.056
    v0_p_cut = v0_p < 1.2*1.056
    chi2_cut = (electron_chi2 < 40) & (positron_chi2 < 40)

    top_track_cluster_dt = top_cluster_time - top_time
    bot_track_cluster_dt = bot_cluster_time - bot_time
    track_cluster_dt_cut = ((np.absolute(top_track_cluster_dt - 43) < 4.5) 
                            & (np.absolute(bot_track_cluster_dt - 43) < 4.5))
    cluster_time_diff_cut = np.absolute(cluster_time_diff) < 2
    
    base_selection = radiative_cut & v0_p_cut & chi2_cut & track_cluster_dt_cut & cluster_time_diff_cut

    # WABS
    l1_cut      = (positron_has_l1 == 1)
    l2_cut      = (positron_has_l2 == 1)
    positron_d0_cut = positron_d0 < 1.1
    asym = (electron_pt - positron_pt)/(electron_pt + positron_pt)
    asym_cut = asym < .47

    wab_cuts = l1_cut & l2_cut & positron_d0_cut & asym_cut

    selection = base_selection & wab_cuts

    with PdfPages("trident_selection.pdf") as pdf :
        
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10)) 
        bins = np.linspace(-10, 10, 201)
        ax.hist(cluster_time_diff, bins, histtype="step", lw=1.5, 
                label="Opp. ECal cluster, trk-cluster match $\chi^2 < 10$, $p(e^-) < 0.75E_{beam}$")
        ax.hist(cluster_time_diff[radiative_cut], bins, histtype="step", lw=1.5, label="Radiative cut")
        ax.hist(cluster_time_diff[radiative_cut & track_cluster_dt_cut],
                bins, histtype="step", lw=1.5, label="abs((Ecal cluster time - track time) - 43) $<$ 4.5")
        ax.hist(cluster_time_diff[radiative_cut & track_cluster_dt_cut & v0_p_cut],
                bins, histtype="step", lw=1.5, label="$p(V_0) < 1.2E_{beam}$")
        ax.hist(cluster_time_diff[radiative_cut  & track_cluster_dt_cut & v0_p_cut & chi2_cut],
                bins, histtype="step", lw=1.5, label="track $\chi^2 < 40$")
        ax.hist(cluster_time_diff[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut],
                bins, histtype="step", lw=1.5,
                label="Top cluster time - Bot cluster time $<$ 2 ns")
        ax.hist(cluster_time_diff[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="l1 and l2 hit")
        ax.hist(cluster_time_diff[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$d0(e^+) < 1.1$")
        ax.hist(cluster_time_diff[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                  & asym_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$p_t(e^-) - p_t(e^+)/p_t(e^-) + p_t(e^+) < .47$")
        ax.set_xlabel("Top cluster time - Bottom cluster time (ns)")
        ax.set_yscale("symlog")
        ax.legend(loc=2, framealpha=0.0, frameon=False)
        pdf.savefig()
        plt.close()


        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10)) 
        bins = np.linspace(0, 1.5, 151)
        ax.hist(v0_p, bins, histtype="step", lw=1.5, 
                label="Opp. ECal cluster, trk-cluster match $\chi^2 < 10$, $p(e^-) < 0.75E_{beam}$")
        ax.hist(v0_p[radiative_cut], bins, histtype="step", lw=1.5, label="Radiative cut")
        ax.hist(v0_p[radiative_cut & track_cluster_dt_cut],
                bins, histtype="step", lw=1.5, label="abs((Ecal cluster time - track time) - 43) $<$ 4.5")
        ax.hist(v0_p[radiative_cut & track_cluster_dt_cut & v0_p_cut],
                bins, histtype="step", lw=1.5, label="$p(V_0) < 1.2E_{beam}$")
        ax.hist(v0_p[radiative_cut  & track_cluster_dt_cut & v0_p_cut & chi2_cut],
                bins, histtype="step", lw=1.5, label="track $\chi^2 < 40$")
        ax.hist(v0_p[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut],
                bins, histtype="step", lw=1.5,
                label="Top cluster time - Bot cluster time $<$ 2 ns")
        ax.hist(v0_p[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="l1 and l2 hit")
        ax.hist(v0_p[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$d0(e^+) < 1.1$")
        ax.hist(v0_p[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                  & asym_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$p_t(e^-) - p_t(e^+)/p_t(e^-) + p_t(e^+) < .47$")
        ax.set_xlabel("$p(V_{0})$ (GeV)")
        ax.legend(loc=2, framealpha=0.0, frameon=False)
        pdf.savefig()
        plt.close()

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10)) 
        bins = np.linspace(0, 1.5, 151)
        ax.hist(electron_p, bins, histtype="step", lw=1.5, 
                label="Opp. ECal cluster, trk-cluster match $\chi^2 < 10$, $p(e^-) < 0.75E_{beam}$")
        ax.hist(electron_p[radiative_cut], bins, histtype="step", lw=1.5, label="Radiative cut")
        ax.hist(electron_p[radiative_cut & track_cluster_dt_cut],
                bins, histtype="step", lw=1.5, label="abs((Ecal cluster time - track time) - 43) $<$ 4.5")
        ax.hist(electron_p[radiative_cut & track_cluster_dt_cut & v0_p_cut],
                bins, histtype="step", lw=1.5, label="$p(V_0) < 1.2E_{beam}$")
        ax.hist(electron_p[radiative_cut  & track_cluster_dt_cut & v0_p_cut & chi2_cut],
                bins, histtype="step", lw=1.5, label="track $\chi^2 < 40$")
        ax.hist(electron_p[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut],
                bins, histtype="step", lw=1.5,
                label="Top cluster time - Bot cluster time $<$ 2 ns")
        ax.hist(electron_p[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="l1 and l2 hit")
        ax.hist(electron_p[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$d0(e^+) < 1.1$")
        ax.hist(electron_p[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                  & asym_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$p_t(e^-) - p_t(e^+)/p_t(e^-) + p_t(e^+) < .47$")
        ax.set_xlabel("$p(e^-)$ (GeV)")
        ax.legend(loc=1, framealpha=0.0, frameon=False)
        pdf.savefig()
        plt.close()

        
        top_track_cluster_dt = np.absolute(top_track_cluster_dt - 43)
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10)) 
        bins = np.linspace(0, 60, 121)
        ax.hist(top_track_cluster_dt, bins, histtype="step", lw=1.5, 
                label="Opp. ECal cluster, trk-cluster match $\chi^2 < 10$, $p(e^-) < 0.75E_{beam}$")
        ax.hist(top_track_cluster_dt[radiative_cut], bins, histtype="step", lw=1.5, label="Radiative cut")
        ax.hist(top_track_cluster_dt[radiative_cut & track_cluster_dt_cut],
                bins, histtype="step", lw=1.5, label="abs((Ecal cluster time - track time) - 43) $<$ 4.5")
        ax.hist(top_track_cluster_dt[radiative_cut & track_cluster_dt_cut & v0_p_cut],
                bins, histtype="step", lw=1.5, label="$p(V_0) < 1.2E_{beam}$")
        ax.hist(top_track_cluster_dt[radiative_cut  & track_cluster_dt_cut & v0_p_cut & chi2_cut],
                bins, histtype="step", lw=1.5, label="track $\chi^2 < 40$")
        ax.hist(top_track_cluster_dt[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut],
                bins, histtype="step", lw=1.5,
                label="Top cluster time - Bot cluster time $<$ 2 ns")
        ax.hist(top_track_cluster_dt[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="l1 and l2 hit")
        ax.hist(top_track_cluster_dt[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$d0(e^+) < 1.1$")
        ax.hist(top_track_cluster_dt[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                  & asym_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$p_t(e^-) - p_t(e^+)/p_t(e^-) + p_t(e^+) < .47$")
        ax.set_xlabel("abs(ECal cluster time - track time - 43) ns")
        ax.legend(loc=1, framealpha=0.0, frameon=False)
        ax.set_yscale("symlog")
        pdf.savefig()
        plt.close()


        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10)) 
        bins = np.linspace(-20, 20, 201)
        ax.hist(positron_d0, bins, histtype="step", lw=1.5, 
                label="Opp. ECal cluster, trk-cluster match $\chi^2 < 10$, $p(e^-) < 0.75E_{beam}$")
        ax.hist(positron_d0[radiative_cut], bins, histtype="step", lw=1.5, label="Radiative cut")
        ax.hist(positron_d0[radiative_cut & track_cluster_dt_cut],
                bins, histtype="step", lw=1.5, label="abs((Ecal cluster time - track time) - 43) $<$ 4.5")
        ax.hist(positron_d0[radiative_cut & track_cluster_dt_cut & v0_p_cut],
                bins, histtype="step", lw=1.5, label="$p(V_0) < 1.2E_{beam}$")
        ax.hist(positron_d0[radiative_cut  & track_cluster_dt_cut & v0_p_cut & chi2_cut],
                bins, histtype="step", lw=1.5, label="track $\chi^2 < 40$")
        ax.hist(positron_d0[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut],
                bins, histtype="step", lw=1.5,
                label="Top cluster time - Bot cluster time $<$ 2 ns")
        ax.hist(positron_d0[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="l1 and l2 hit")
        ax.hist(positron_d0[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$d0(e^+) < 1.1$")
        ax.hist(positron_d0[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                  & asym_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$p_t(e^-) - p_t(e^+)/p_t(e^-) + p_t(e^+) < .47$")
        ax.set_xlabel("Positron D0 (mm)")
        ax.legend(loc=1, framealpha=0.0, frameon=False)
        ax.set_yscale("symlog")
        pdf.savefig()
        plt.close()

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10)) 
        bins = np.linspace(-1, 1, 201)
        ax.hist(asym, bins, histtype="step", lw=1.5, 
                label="Opp. ECal cluster, trk-cluster match $\chi^2 < 10$, $p(e^-) < 0.75E_{beam}$")
        ax.hist(asym[radiative_cut], bins, histtype="step", lw=1.5, label="Radiative cut")
        ax.hist(asym[radiative_cut & track_cluster_dt_cut],
                bins, histtype="step", lw=1.5, label="abs((Ecal cluster time - track time) - 43) $<$ 4.5")
        ax.hist(asym[radiative_cut & track_cluster_dt_cut & v0_p_cut],
                bins, histtype="step", lw=1.5, label="$p(V_0) < 1.2E_{beam}$")
        ax.hist(asym[radiative_cut  & track_cluster_dt_cut & v0_p_cut & chi2_cut],
                bins, histtype="step", lw=1.5, label="track $\chi^2 < 40$")
        ax.hist(asym[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut],
                bins, histtype="step", lw=1.5,
                label="Top cluster time - Bot cluster time $<$ 2 ns")
        ax.hist(asym[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="l1 and l2 hit")
        ax.hist(asym[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$d0(e^+) < 1.1$")
        ax.hist(asym[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                  & asym_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$p_t(e^-) - p_t(e^+)/p_t(e^-) + p_t(e^+) < .47$")
        ax.set_xlabel("$p_t(e^-) - p_t(e^+)/p_t(e^-) + p_t(e^+)$")
        ax.legend(loc=1, framealpha=0.0, frameon=False)
        ax.set_yscale("symlog")
        pdf.savefig()
        plt.close()


        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10)) 
        bins = np.linspace(0, 100, 201)
        ax.hist(electron_chi2, bins, histtype="step", lw=1.5, 
                label="Opp. ECal cluster, trk-cluster match $\chi^2 < 10$, $p(e^-) < 0.75E_{beam}$")
        ax.hist(electron_chi2[radiative_cut], bins, histtype="step", lw=1.5, label="Radiative cut")
        ax.hist(electron_chi2[radiative_cut & track_cluster_dt_cut],
                bins, histtype="step", lw=1.5, label="abs((Ecal cluster time - track time) - 43) $<$ 4.5")
        ax.hist(electron_chi2[radiative_cut & track_cluster_dt_cut & v0_p_cut],
                bins, histtype="step", lw=1.5, label="$p(V_0) < 1.2E_{beam}$")
        ax.hist(electron_chi2[radiative_cut  & track_cluster_dt_cut & v0_p_cut & chi2_cut],
                bins, histtype="step", lw=1.5, label="track $\chi^2 < 40$")
        ax.hist(electron_chi2[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut],
                bins, histtype="step", lw=1.5,
                label="Top cluster time - Bot cluster time $<$ 2 ns")
        ax.hist(electron_chi2[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="l1 and l2 hit")
        ax.hist(electron_chi2[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$d0(e^+) < 1.1$")
        ax.hist(electron_chi2[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                  & asym_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$p_t(e^-) - p_t(e^+)/p_t(e^-) + p_t(e^+) < .47$")
        ax.set_xlabel("Track $\chi^2$")
        ax.legend(loc=1, framealpha=0.0, frameon=False)
        ax.set_yscale("symlog")
        pdf.savefig()
        plt.close()


        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10)) 
        bins = np.linspace(0, .1, 2000)
        ax.hist(mass, bins, histtype="step", lw=1.5, 
                label="Opp. ECal cluster, trk-cluster match $\chi^2 < 10$, $p(e^-) < 0.75E_{beam}$")
        ax.hist(mass[radiative_cut], bins, histtype="step", lw=1.5, label="Radiative cut")
        ax.hist(mass[radiative_cut & track_cluster_dt_cut],
                bins, histtype="step", lw=1.5, label="abs((Ecal cluster time - track time) - 43) $<$ 4.5")
        ax.hist(mass[radiative_cut & track_cluster_dt_cut & v0_p_cut],
                bins, histtype="step", lw=1.5, label="$p(V_0) < 1.2E_{beam}$")
        ax.hist(mass[radiative_cut  & track_cluster_dt_cut & v0_p_cut & chi2_cut],
                bins, histtype="step", lw=1.5, label="track $\chi^2 < 40$")
        ax.hist(mass[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut],
                bins, histtype="step", lw=1.5,
                label="Top cluster time - Bot cluster time $<$ 2 ns")
        ax.hist(mass[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="l1 and l2 hit")
        ax.hist(mass[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$d0(e^+) < 1.1$")
        ax.hist(mass[radiative_cut  
                                  & track_cluster_dt_cut 
                                  & v0_p_cut 
                                  & chi2_cut 
                                  & cluster_time_diff_cut
                                  & l1_cut & l2_cut
                                  & positron_d0_cut
                                  & asym_cut
                                 ],
                bins, histtype="step", lw=1.5,
                label="$p_t(e^-) - p_t(e^+)/p_t(e^-) + p_t(e^+) < .47$")
        ax.set_xlabel("Invariant mass (GeV)")
        ax.legend(loc=1, framealpha=0.0, frameon=False)
        ax.set_yscale("symlog")
        ax.set_xlim(0, .1)
        pdf.savefig()
        plt.close()
           
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10)) 
        bins = np.linspace(0, .1, 2000)
        ax.hist(mass[selection], bins, histtype="step", lw=1.5)
        ax.set_xlabel("Invariant mass (GeV)")
        ax.set_yscale("symlog")
        ax.legend(loc=1)
        pdf.savefig()
        plt.close()

        file = r.TFile("invariant_mass_final_selection.root", "recreate")
        mass_histo = r.TH1F("invariant_mass", "invariant_mass", 2000, 0., 0.1)
        #mass_histo = r.TH1F("invariant_mass", "invariant_mass", 50, .0, 0.1)
        for value in np.nditer(mass[selection]) : 
            mass_histo.Fill(value)
        
        mass_histo.Write()
        file.Close()



if __name__ == "__main__" : 
    main()

'''
class TridentSelection : 

    def process( :

        chi2_cut      = (electron_chi2 < 40) & (positron_chi2 < 40)
        fee_cut       = electron_p < 0.8
        p_sum_cut     = v0_p < 1.2
        radiative_cut = v0_p > 0.8
        selection     = chi2_cut & fee_cut & p_sum_cut & radiative_cut

        with PdfPages("trident_selection.pdf") as pdf :

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))  
            bins = np.linspace(0, 100, 101)
            ax.hist(electron_chi2, bins, histtype="step")
            ax.hist(electron_chi2[chi2_cut], bins, histtype="step")
            ax.set_xlabel("$e^-$ $\chi^2$")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))  
            bins = np.linspace(0, 100, 101)
            ax.hist(positron_chi2, bins, histtype="step")
            ax.hist(positron_chi2[chi2_cut], bins, histtype="step")
            ax.set_xlabel("$e^+$ $\chi^2$")
            pdf.savefig()
            plt.close()

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))   

            bins = np.linspace(0, 1.5, 151)
            ax.hist(electron_p, bins, histtype="step", lw=1.5, label="All")
            ax.hist(electron_p[fee_cut], bins, histtype="step", lw=1.5, label="$e(p) < 0.8")
            pdf.savefig()
            plt.close()
            '''

