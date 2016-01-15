
#ifndef __BUMP_HUNTER_H__
#define __BUMP_HUNTER_H__

#include <vector>
#include <map>

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

        std::map <std::string, RooRealVar*> variable_map; 

        RooAddPdf* model;  
        RooGaussian* signal;
        RooChebychev* bkg; 
        RooArgList arg_list;

        double window_size; 
};

#endif // __BUMP_HUNTER_H__
