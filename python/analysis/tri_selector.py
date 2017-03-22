#!/usr/bin/python

from __future__ import division

import argparse
import collections
import sys
import math
import numpy as np
import Plotter
import root_numpy as rnp
import ROOT as r

def main() :

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-f", "--file_list", help="List of ROOT files to process.")
    parser.add_argument("-l", "--lumi",      help="Luminosity")
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

    apply_tri_selection(rec, args.lumi)

def apply_tri_selection(rec, lumi):

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

    top_track_cluster_dt = top_cluster_time - top_time
    abs_top_track_cluster_dt = np.absolute(top_track_cluster_dt - 43)
    bot_track_cluster_dt = bot_cluster_time - bot_time
    abs_bot_track_cluster_dt = np.absolute(bot_track_cluster_dt - 43)
    track_cluster_dt_cut = ((abs_top_track_cluster_dt < 4.5) 
                            & (abs_bot_track_cluster_dt < 4.5))

    asym = (electron_pt - positron_pt)/(electron_pt + positron_pt)
    #
    # Define cuts
    #
    cuts = collections.OrderedDict()

    # Base cuts used to reduce accidentals
    cuts['Radiative cut'] = v0_p > 0.8*1.056 # GeV
    cuts['abs(Ecal clust time - trk time) - 43 ns < 4.5'] = track_cluster_dt_cut
    cuts['$p(V_0) < 1.2 E_{beam}$'] = v0_p < 1.2*1.056 # GeV
    cuts['trk $\chi^2$ < 40'] = (electron_chi2 < 40) & (positron_chi2 < 40)
    cuts['Ecal clust pair dt < 2 ns'] = np.absolute(cluster_time_diff) < 2
    cuts['l1 & l2 hit'] = (positron_has_l1 == 1) & (positron_has_l2 == 1)
    cuts['$d_{0}(e^+) < 1.1$'] = positron_d0 < 1.1
    cuts['$p_t(e^-) - p_t(e^+)/p_t(e^-) + p_t(e^+)$'] = asym < .47
    
    labels = ['Opp. Ecal clusters, trk-cluster match $\chi^2 < 10$, $p(e^-)<0.75E_{beam}$']
    clust_dt_arr = [cluster_time_diff]
    v0_p_arr = [v0_p]
    electron_p_arr = [electron_p]
    positron_d0_arr = [positron_d0]
    electron_chi2_arr = [electron_chi2]
    asym_arr = [asym]
    abs_top_cluster_dt_arr = [abs_top_track_cluster_dt]
    
    cut = np.ones(len(v0_p), dtype=bool)
    for key, value in cuts.iteritems():
        cut = cut & value
        clust_dt_arr.append(cluster_time_diff[cut])
        v0_p_arr.append(v0_p[cut])
        electron_p_arr.append(electron_p[cut])
        positron_d0_arr.append(positron_d0[cut])
        electron_chi2_arr.append(electron_chi2[cut])
        abs_top_cluster_dt_arr.append(abs_top_track_cluster_dt[cut])
        asym_arr.append(asym[cut])
        labels.append(key)

    plt = Plotter.Plotter('trident_selection.pdf')
    
    plt.plot_hists(clust_dt_arr, 
                   np.linspace(-10, 10, 201),
                   labels=labels,
                   x_label='Top cluster time - Bottom cluster time (ns)',
                   label_loc=2, 
                   norm=True,
                   ylog=True,
                   root=True)

    plt.plot_hists(v0_p_arr, 
                   np.linspace(0, 1.5, 151),
                   labels=labels,
                   label_loc=2,
                   x_label='$V_{0}(p)$ (GeV)',
                   ylog=True)
    
    plt.plot_hists(electron_p_arr, 
                   np.linspace(0, 1.5, 151), 
                   labels=labels, 
                   x_label='$p(e^-)$ (GeV)',
                   norm=True,
                   label_loc=4, 
                   ylog=True)

    plt.plot_hists(positron_d0_arr, 
                   np.linspace(-20, 20, 201),
                   labels=labels, 
                   x_label='$d_{0}(e^{+})$', 
                   norm=True,
                   label_loc=2, 
                   ylog=True)

    plt.plot_hists(electron_chi2_arr, 
                   np.linspace(0, 100, 201),
                   labels=labels, 
                   x_label='Track $\chi^2$', 
                   norm=True,
                   label_loc=1, 
                   ylog=True)

    plt.plot_hists(asym_arr, 
                   np.linspace(-1, 1, 201),
                   labels=labels, 
                   x_label='$p_t(e^-) - p_t(e^+)/p_t(e^-) + p_t(e^+)$',
                   norm=True,
                   label_loc=2, 
                   ylog=True)

    plt.plot_hists(abs_top_cluster_dt_arr,
                   np.linspace(0, 60, 121),
                   labels=labels,
                   x_label='abs(ECal cluster time - track time - 43) ns',
                   norm=True,
                   label_loc=1,
                   ylog=True)

    plt.close()

    file = r.TFile("invariant_mass.root", "recreate")

    mass_histo = r.TH1F("invariant_mass", "invariant_mass", 2000, 0., 0.1)
    #mass_histo = r.TH1F("invariant_mass", "invariant_mass", 50, 0., 0.1)
    mass_histo.GetXaxis().SetTitle("m(e^+e^-) (GeV)")
    mass_histo.GetYaxis().SetTitle("#sigma(#mub)")
    bin_width = mass_histo.GetXaxis().GetBinWidth(1)
        
    weights = np.empty(len(mass[cut]))
    if lumi: weights.fill(1/(bin_width*float(lumi)))
    else: weights.fill(1)
        #rnp.fill_hist(mass_histo, mass[selection], weights=weights)
    rnp.fill_hist(mass_histo, mass[cut], weights=weights)    
    mass_histo.Write()
    file.Close()


'''
        
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

        #mass_histo = r.TH1F("invariant_mass", "invariant_mass", 2000, 0., 0.1)
        mass_histo = r.TH1F("invariant_mass", "invariant_mass", 50, 0., 0.1)
        mass_histo.GetXaxis().SetTitle("m(e^+e^-) (GeV)")
        mass_histo.GetYaxis().SetTitle("#sigma(#mub)")
        bin_width = mass_histo.GetXaxis().GetBinWidth(1)
        
        weights = np.empty(len(mass[cut]))
        if lumi: weights.fill(1/(bin_width*float(lumi)))
        else: weights.fill(1)
        #rnp.fill_hist(mass_histo, mass[selection], weights=weights)
        rnp.fill_hist(mass_histo, mass[cut], weights=weights)
        
        mass_histo.Write()
        file.Close()
'''


if __name__ == "__main__" : 
    main()
