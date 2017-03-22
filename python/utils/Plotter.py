
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import ROOT as r
import root_numpy as rnp

from itertools import izip
from matplotlib.backends.backend_pdf import PdfPages

class Plotter(object):

    def __init__(self, file_path): 
        
        plt.style.use('bmh')
        matplotlib.rcParams.update({'font.size': 12})
        matplotlib.rcParams['axes.facecolor'] = 'white'
        matplotlib.rcParams['legend.numpoints'] = 1
        
        self.pdf = PdfPages(file_path + '.pdf')
        print '[ Plotter ] Saving plots to %s' % (file_path + '.pdf')

        self.rfile = r.TFile(file_path + '.root', 'recreate')
        print '[ Plotter ] Writing histograms to %s' % (file_path + '.root')

    def plot_hist(self, values, bins, **params):
       
        name = None
        if 'name' in params: 
            name = params['name']

        if ('root' in params): 
            histo = r.TH1F(name, name, len(bins), np.amin(bins), np.amax(bins))
            rnp.fill_hist(histo, values)
            histo.Scale(1/histo.Integral(), 'width')
            histo.Write()
            
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
        
        label=None
        if 'label' in params:
            label=params['label']
            ax.legend()
       
        norm = False
        if 'norm' in params: 
            norm = True

        if 'ylog' in params:
            if norm: ax.set_yscale('log')
            else: ax.set_yscale('symlog')

        if 'x_label' in params:
            ax.set_xlabel(params['x_label'])


        ax.hist(values, bins, histtype='step', lw=1.5, label=label, normed=norm)

        self.pdf.savefig(bbox_inches='tight')
        plt.close()

    def plot_hists(self, values, bins, **params):
    
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))

        label=None

        norm = False
        if 'norm' in params: 
            norm = True

        if 'ylog' in params:
            if norm: ax.set_yscale('log')
            else: ax.set_yscale('symlog')

        if 'x_label' in params:
            ax.set_xlabel(params['x_label'])

        labels = None
        if 'labels' in params:
            labels = params['labels']
 
        label_loc=0
        if 'label_loc' in params: 
            label_loc=params['label_loc']

        for x_arr, label in izip(values, labels):
            ax.hist(x_arr, bins, histtype='step', lw=1.5, normed=norm, label=label)

        if labels: ax.legend(framealpha=0.0, frameon=False, loc=label_loc)

        self.pdf.savefig(bbox_inches='tight')
        plt.close()

        if ('root' in params): 
            for x_arr, label in izip(values, labels):
                histo = r.TH1F(label, label, len(bins), np.amin(bins), np.amax(bins))
                rnp.fill_hist(histo, x_arr)
                histo.Scale(1/histo.Integral(), 'width')
                histo.Write()


    def plot_graph(self, x, y, x_err, y_err, **params):
       
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
        
        label=None
        if 'label' in params:
            label=params['label']
            ax.legend()
       
        if 'x_label' in params:
            ax.set_xlabel(params['x_label'])

        if 'y_label' in params:
            ax.set_ylabel(params['y_label'])

        if 'ylim' in params: 
            ax.set_ylim(params['ylim'])

        if 'xlog' in params:
            ax.set_xscale('symlog')

        ax.errorbar(x, y, x_err, y_err, markersize=10, marker='o', 
                    linestyle='-', fmt='', label=label)


        self.pdf.savefig(bbox_inches='tight')
        plt.close()


    def plot_graphs(self, x, y, x_err, y_err, **params):
       
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
        
        labels=None
        if 'labels' in params:
            labels=params['labels']

        label_loc=0
        if 'label_loc' in params: 
            label_loc=params['label_loc']
       
        if 'x_label' in params:
            ax.set_xlabel(params['x_label'])

        if 'y_label' in params:
            ax.set_ylabel(params['y_label'])

        if 'ylim' in params: 
            ax.set_ylim(params['ylim'])
        
        if 'xlog' in params:
            ax.set_xscale('symlog')

        for index in xrange(0, len(x)):
            ax.errorbar(x[index], y[index], 0, 0, 
                        markersize=6, marker='o', 
                        linestyle='-', fmt='', label=labels[index])

        if labels: ax.legend(loc=label_loc)
        
        self.pdf.savefig(bbox_inches='tight')
        plt.close()


    def close(self):
        self.pdf.close()
        self.rfile.Close()



        

