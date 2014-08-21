
#
#	@author: Omar Moreno <omoreno1@ucsc.edu>
#
#

from ROOT import RooGaussian

class RooFitResult(object):

	def __init__(self, pdf, fit_params, fit_param_errors):
		self.pdf = 	pdf
		self.fit_params = fit_params
		self.fit_param_errors = fit_param_errors

	def get_pdf(self):
		return self.pdf

class RooGaussFitResult(RooFitResult):

	mean_index = 0
	sigma_index = 1

	def get_mean(self):
		return self.fit_params[RooGaussFitResult.mean_index]

	def get_sigma(self):
		return self.fit_params[RooGaussFitResult.sigma_index]

	def get_mean_error(self):
		return self.fit_param_errors[RooGaussFitResult.mean_index]

	def get_sigma_error(self):
		return self.fit_param_errors[RooGaussFitResult.sigma_index]


