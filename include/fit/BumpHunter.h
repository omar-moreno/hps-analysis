
#ifndef __BUMP_HUNTER_H__
#define __BUMP_HUNTER_H__

#include <vector>

#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooArgList.h>
#include <RooAddPdf.h>

class BumpHunter {

    public:

        /** Default Constructor */
        BumpHunter(int poly_order);

        /** Destructor */
        ~BumpHunter();

    private: 

        RooRealVar* invariant_mass;
        RooRealVar* ap_mass_mean;
        RooRealVar* ap_mass_sigma;
        RooRealVar* n_sig;
        RooRealVar* n_bkg;
        RooAddPdf* model;  
        RooGaussian* signal;
        RooChebychev* bkg; 
        RooArgList arg_list;

        std::vector<RooRealVar*> t;

        double window_size; 
};

#endif // __BUMP_HUNTER_H__
