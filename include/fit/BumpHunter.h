/**
 * @file BumpHunter.h
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date January 14, 2015
 */

#ifndef __BUMP_HUNTER_H__
#define __BUMP_HUNTER_H__

//----------------//   
//   C++ StdLib   //
//----------------//   
#include <vector>
#include <map>
#include <fstream>

//----------//
//   ROOT   //
//----------//   
#include <TH1.h>

//------------//
//   RooFit   //
//------------//
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
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
         * Fit the given histogram in the range window_start, window_end.
         *
         * @param histogram The histogram to fit.
         * @param window_start
         * @param window_end
         * @param window_step 
         */
        std::map<double, RooFitResult*> fit(TH1* histogram, double window_start, double window_end, double window_step);

        /**
         * 
         */
        void setWindowSize(double window_size) { this->window_size = window_size; }; 

        /** Fit using a background only model. */
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
   
        /**
         * Reset the fit parameters to their initial values.
         *
         * @param initial_params A list containing the fit parameters.
         */ 
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

        /** Output file stream */
        std::ofstream* ofs;

        double window_size;

        bool bkg_only;  
};

#endif // __BUMP_HUNTER_H__
