
#
#   @author: Omar Moreno <omoreno1@ucsc.edu>
#
#

import ROOT as r

from ROOT import RooRealVar, RooGaussian
from RooFitResult import RooGaussFitResult
from ROOT import RooAbsPdf

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

    def fit_to_gaus_exp(self):
        histogram_min = self.histogram.GetXaxis().GetXmin()
        histogram_max = self.histogram.GetXaxis().GetXmax()
        x = r.RooRealVar("x", "x", histogram_min, histogram_max)
        plot = x.frame()

        arg_list = r.RooArgList(x)
        histogram_data = r.RooDataHist("histogram_data", "histogram_data", arg_list, self.histogram)
        mean = histogram_data.mean(x)
        sigma = histogram_data.sigma(x)
       
        
        mean_var = RooRealVar("mean_var", "Gaussian - Mean", mean, mean - 10*sigma, mean + 10*sigma)
        sigma_var = RooRealVar("sigma_var", "Gaussian - Sigma", sigma, 0, 10*sigma)
        gaussian = RooGaussian("gaussian", "Gaussian", x, mean_var, sigma_var)
        
        a_var = RooRealVar("rlife", "rlife", 0, -10000, 0)
        exponential = r.RooExponential("expo", "expo", x, a_var)
        
        print "Mean: " + str(mean)
        print "Sigma: " + str(sigma)
        pdf = r.RooNumConvPdf("gaus_x_exp", "gaus_x_exp", x, gaussian, exponential)
        pdf.setConvolutionWindow(mean_var, sigma_var, 5)

        #pdf = r.RooGExpModel("pdf", "pdf", x, sigma_var, rlife_var)
        pdf.fitTo(histogram_data)
        histogram_data.plotOn(plot)
        gaussian.plotOn(plot, r.RooFit.LineStyle(r.kDashed), r.RooFit.LineColor(r.kRed))
        exponential.plotOn(plot, r.RooFit.LineColor(r.kGreen))
        pdf.plotOn(plot)

        plot.Draw()
        


#   def fit_to_gaussian(self, data, variable, plotter):
#
#       data.plotOn(plotter)
#
#       mean = data.mean(variable)
#       sigma = data.sigma(variable)
        
#       mean_var = RooRealVar("mean_var", "Gaussian - Mean", mean, mean - 4*sigma, mean + 4*sigma)
#       sigma_var = RooRealVar("sigma_var", "Gaussian - Sigma", sigma, 0, 10*sigma)
        
#       gaussian = RooGaussian("gaussian", "Gaussian", variable, mean_var, sigma_var)
#       fit_result = gaussian.fitTo(data)
#       gaussian.plotOn(plotter)

#       fit_params = [0]*2
#       fit_params[0] = mean_var.getValV()
#       fit_params[1] = sigma_var.getValV()
#       fit_param_errors = [0]*2
#       fit_param_errors[0] = mean_var.getError()
#       fit_param_errors[1] = sigma_var.getError()

#       return RooGaussFitResult(gaussian, fit_params, fit_param_errors) 

