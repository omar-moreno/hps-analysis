
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
#include <RooMinuit.h>
#include <RooFitResult.h>

class BumpHunter {

    public:

        /** Default Constructor */
        BumpHunter(int poly_order);

        /** Destructor */
        ~BumpHunter();

        /**
         *
         */
        std::map<double, RooFitResult*> fit(TH1* histogram, double window_start, double window_end, double window_step);

        void setWindowSize(double window_size) { this->window_size = window_size; }; 

        void fitBkgOnly();

    private: 

        /**
         * Get the HPS mass resolution at the given mass.  The functional form 
         * of the mass resolution was determined using MC.
         *
         * @param mass The mass of interest.
         * @return The mass resolution at the given mass.
         */
        inline double getMassResolution(double mass) { 
            return -6.166*mass*mass*mass + 0.9069*mass*mass -0.00297*mass + 0.000579; 
        };
    
        void resetParameters(RooArgList initial_params); 

        std::map <std::string, RooRealVar*> variable_map; 

        /** Signal + bkg model */
        RooAddPdf* comp_model;  

        /** Bkg only model */
        RooAddPdf* bkg_model;

        /** */
        RooAddPdf* model; 

        /** Signal PDF */ 
        RooGaussian* signal;

        /** Bkg PDF */
        RooChebychev* bkg;

        /** */ 
        RooArgList arg_list;

        double window_size;

        bool bkg_only;  
};

#endif // __BUMP_HUNTER_H__
