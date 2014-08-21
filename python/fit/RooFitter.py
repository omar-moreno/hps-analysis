
#
#	@author: Omar Moreno <omoreno1@ucsc.edu>
#
#

from ROOT import RooRealVar, RooGaussian
from RooFitResult import RooGaussFitResult
from ROOT import RooAbsPdf

class RooFitter(object):

	def fit_to_gaussian(self, data, variable, plotter):

		data.plotOn(plotter)

		mean = data.mean(variable)
		sigma = data.sigma(variable)
		
		mean_var = RooRealVar("mean_var", "Gaussian - Mean", mean, mean - 4*sigma, mean + 4*sigma)
		sigma_var = RooRealVar("sigma_var", "Gaussian - Sigma", sigma, 0, 10*sigma)
		
		gaussian = RooGaussian("gaussian", "Gaussian", variable, mean_var, sigma_var)
		fit_result = gaussian.fitTo(data)
		gaussian.plotOn(plotter)

		fit_params = [0]*2
		fit_params[0] = mean_var.getValV()
		fit_params[1] = sigma_var.getValV()
		fit_param_errors = [0]*2
		fit_param_errors[0] = mean_var.getError()
		fit_param_errors[1] = sigma_var.getError()

		return RooGaussFitResult(gaussian, fit_params, fit_param_errors) 


