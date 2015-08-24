/**
 * @file RooFitter.h
 * @brief A common set of fits done using RooFit
 * @author Omar Moreno <omoreno1@ucsc.edu> \n
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date August 24, 2015
 */

#ifndef __ROO_FITTER_H__
#define __ROO_FITTER_H__

#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooAddPdf.h>

#include <TH1.h>

namespace RooFitter { 

    RooPlot* fitToGaussian(TH1* histogram, RooRealVar var); 

    RooPlot* fitToDoubleGaussian(TH1* histogram, RooRealVar var); 
}

#endif // __ROO_FITTER_H__
