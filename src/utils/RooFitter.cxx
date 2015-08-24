/**
 * @file RooFitter.cxx
 * @brief A common set of fits done using RooFit
 * @author Omar Moreno <omoreno1@ucsc.edu> \n
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date August 24, 2015
 */

#include <RooFitter.h>

using namespace RooFit; 

RooPlot* RooFitter::fitToGaussian(TH1* histogram, RooRealVar var) { 

    RooDataHist histogram_data("data", histogram->GetName(), var, histogram); 

    double mean = histogram->GetMean(); 
    double rms = histogram->GetRMS(); 

    RooRealVar gauss_mean("gauss_mean", "mean", mean, mean - 4*rms, mean + 4*rms);  
    RooRealVar gauss_sigma("gauss_sigma", "#sigma", rms, rms - 3*rms, rms + 3*rms); 
    RooGaussian gaussian("gaussian", "gaussian", var, gauss_mean, gauss_sigma); 

    RooFitResult* result = gaussian.fitTo(histogram_data, Range(mean - 2.*rms, mean + 2.*rms), Save());

    RooPlot* plot = var.frame(); 
    histogram_data.plotOn(plot); 
    gaussian.plotOn(plot, Range(0, 5.0));
    gaussian.paramOn(plot);   

    return plot; 
}

RooPlot* RooFitter::fitToDoubleGaussian(TH1* histogram, RooRealVar var) { 
    
    RooDataHist histogram_data("data", histogram->GetName(), var, histogram); 

    double mean = histogram->GetMean(); 
    double rms = histogram->GetRMS(); 

    RooRealVar sig_gauss_mean("sig_gauss_mean", "mean", mean, mean - 3*rms, mean + 3*rms);  
    RooRealVar sig_gauss_sigma("sig_gauss_sigma", "#sigma", rms, rms - 3*rms, rms + 3*rms); 
    RooGaussian sig_gaussian("sig_gaussian", "gaussian", var, sig_gauss_mean, sig_gauss_sigma); 

    RooRealVar bkg_gauss_mean("bkg_gauss_mean", "mean", mean, mean - 10*rms, mean + 10*rms);  
    RooRealVar bkg_gauss_sigma("bkg_gauss_sigma", "#sigma", rms, rms - 10*rms, rms + 10*rms); 
    RooGaussian bkg_gaussian("bkg_gaussian", "gaussian", var, bkg_gauss_mean, bkg_gauss_sigma); 

    RooRealVar n_sig("n_sig", "signal events", 5000, 0, 100000);
    RooRealVar n_bkg("n_bkg", "background events", 500, 0, 10000); 
    RooAddPdf sum("sum", "gauss+gauss", RooArgList(sig_gaussian, bkg_gaussian), RooArgList(n_sig, n_bkg)); 
    
    RooFitResult* result = sum.fitTo(histogram_data, Save()); 
    
    RooPlot* plot = var.frame(); 
    histogram_data.plotOn(plot); 
    sig_gaussian.plotOn(plot, LineStyle(kDashed));
    bkg_gaussian.plotOn(plot, LineStyle(kDashed), LineColor(kRed)); 
    sum.plotOn(plot);
    sum.paramOn(plot); 

    return plot;  
}
