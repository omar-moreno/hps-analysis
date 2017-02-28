
import matplotlib
import matplotlib.pyplot as plt

from itertools import izip
from matplotlib.backends.backend_pdf import PdfPages

class Plotter(object):

    def __init__(self, file_path): 
        
        plt.style.use('bmh')
        matplotlib.rcParams.update({'font.size': 8})
        matplotlib.rcParams['axes.facecolor'] = 'white'
        matplotlib.rcParams['legend.numpoints'] = 1
        
        self.pdf = PdfPages(file_path)

    def plot_hist(self, values, bins, **params):
        
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

        label_loc = 1
        if 'label_loc' in params: 
            label_loc = params['label_loc']
 
        for x_arr, label in izip(values, labels):
            ax.hist(x_arr, bins, histtype='step', lw=1.5, normed=norm, label=label)

        if labels: ax.legend(loc=int(label_loc), framealpha=0.0, frameon=False)

        self.pdf.savefig(bbox_inches='tight')
        plt.close()

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
       
        if 'x_label' in params:
            ax.set_xlabel(params['x_label'])

        if 'y_label' in params:
            ax.set_ylabel(params['y_label'])

        if 'ylim' in params: 
            ax.set_ylim(params['ylim'])
        
        if 'xlog' in params:
            ax.set_xscale('symlog')

        for index in xrange(0, len(x)):
            ax.errorbar(x[index], y[index], x_err, y_err, 
                        markersize=6, marker='o', 
                        linestyle='-', fmt='', label=labels[index])

        if labels: ax.legend()
        
        self.pdf.savefig(bbox_inches='tight')
        plt.close()


    def close(self):
        self.pdf.close()



        

