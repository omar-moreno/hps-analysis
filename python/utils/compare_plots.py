
"""
"""

__author__ = "Omar Moreno"
__email__  = "omoreno@slac.stanford.edu"
__date__   = "October 05, 2016"

import argparse
import collections
import matplotlib
import matplotlib.pyplot as plt
import ROOT as r
import root_numpy as rnp
import sys
import yaml

from matplotlib.backends.backend_pdf import PdfPages
from rootpy.io import File
import rootpy.plotting.root2matplotlib as rplt

def parse_config(config_file) : 
    print "Loading configuration from " + str(config_file)
    config = open(config_file, 'r')
    return yaml.load(config)

def root_histo_to_numpy(histogram, normalize) : 
    bin_centers = []
    bin_contents = []
    bin_errors = []

    if normalize :
        histogram.Sumw2()
        histogram.Scale(1/histogram.Integral())

    for bin_index in xrange(0, histogram.GetNbinsX()+1) :
        bin_centers.append(histogram.GetBinCenter(bin_index))
        bin_contents.append(histogram.GetBinContent(bin_index))
        bin_errors.append(histogram.GetBinError(bin_index))
    
    return bin_centers, bin_contents, bin_errors

def main() : 

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-c", "--config", 
                        help="Configuration file.")
    args = parser.parse_args()

    if args.config is None : 
        print 'A configuration file needs to be specified!'
        sys.exit(2)

    plt.style.use('bmh')
    matplotlib.rcParams.update({'font.size': 12})
    matplotlib.rcParams['axes.facecolor'] = 'white'
    matplotlib.rcParams['legend.numpoints'] = 1

    config = parse_config(args.config)

    root_files = []
    for root_file in config["Files"] :
        root_files.append(File.Open(root_file)) 
    
    print root_files

    histograms = collections.OrderedDict()
    for root_file in root_files :
        for histogram in root_file.objects(r.TH1) : 
            if histogram.GetName() not in histograms : 
                histograms[histogram.GetName()] = []

            #print "Histogram Name: " + str(histogram.GetName())
            histograms[histogram.GetName()].append(histogram)
    
    with PdfPages("comparisons.pdf") as pdf : 

        for key, histogram_list in histograms.iteritems() : 
            
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 7.5))
            
            bin_centers_data, bin_contents_data, bin_errors_data = root_histo_to_numpy(histogram_list[0], True)
            bin_centers_mc, bin_contents_mc, bin_errors_mc = root_histo_to_numpy(histogram_list[1], True)
            
            #print "Data: " + str(bin_centers_data)
            #print "MC: " + str(bin_centers_mc)

            plt.errorbar(bin_centers_mc, bin_contents_mc, ls='steps-mid', label="MC")
            plt.errorbar(bin_centers_data, bin_contents_data, yerr=bin_errors_data, xerr=0.5,
                    fmt='o', markersize=3, capsize=0, ls='none', label="Run 5772")
            plt.xlim(-5, 40)
            plt.xlabel("Raw Hit Multiplicity")
            plt.title(key)
            plt.legend()
            
            pdf.savefig(bbox_inches='tight')
            plt.close()
    
if __name__ == "__main__":
    main()
