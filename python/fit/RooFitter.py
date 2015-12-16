
#
#   @author: Omar Moreno <omoreno1@ucsc.edu>
#
#

import ROOT as r

class RooFitter(object):

    def __init__(self, root_file_name, histogram_name) : 
        self.root_file = r.TFile(root_file_name)
        self.histogram = self.root_file.Get(histogram_name)

    def fit_to_gaussian(self) :
        histogram_min = self.histogram.GetXaxis().GetXmin()
        histogram_max = self.histogram.GetXaxis().GetXmax()
        x = r.RooRealVar("x", "x", histogram_min, histogram_max)
        plot = x.frame()
        arg_list = r.RooArgList(x)
        histogram_data = r.RooDataHist("histogram_data", "histogram_data", arg_list, self.histogram)
        mean = histogram_data.mean(x)
        sigma = histogram_data.sigma(x)
        print "Histogram mean: " + str(mean)
        print "Histogram sigma: " + str(sigma)

        mean_var = RooRealVar("mean_var", "Gaussian - Mean", mean, mean - 4*sigma, mean + 4*sigma)
        sigma_var = RooRealVar("sigma_var", "Gaussian - Sigma", sigma, 0, 10*sigma)
        
        gaussian = RooGaussian("gaussian", "Gaussian", x, mean_var, sigma_var)
        fit_result = gaussian.fitTo(histogram_data)
        histogram_data.plotOn(plot)
        gaussian.plotOn(plot)

        plot.Draw()

    def fit_to_crystal_ball(self) : 

        histogram_min = self.histogram.GetXaxis().GetXmin()
        histogram_max = self.histogram.GetXaxis().GetXmax()

        x = r.RooRealVar("x", "x", histogram_min, histogram_max)
        plot = x.frame()

        arg_list = r.RooArgList(x)
        histogram_data = r.RooDataHist("histogram_data", "histogram_data", arg_list, self.histogram)

        mean = histogram_data.mean(x)
        sigma = histogram_data.sigma(x)
        
        mean_var = r.RooRealVar("#mu", "mean", mean, mean - 10*sigma, mean + 10*sigma)
        sigma_var = r.RooRealVar("#sigma", "mean", sigma, 0, 10*sigma)
        alpha_var = r.RooRealVar("#alpha", "alpha", -2, -100, -1)
        n_var = r.RooRealVar("n", "n", 0, 0, 20)

        cb = r.RooCBShape("cb", "cb", x, mean_var, sigma_var, alpha_var, n_var)

        cb.fitTo(histogram_data)

        histogram_data.plotOn(plot)
        cb.plotOn(plot)
        cb.paramOn(plot)

        plot.Draw()
