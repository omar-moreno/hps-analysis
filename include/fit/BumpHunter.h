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
#include <cstdio> 
#include <vector>
#include <map>
#include <fstream>
#include <cmath>
#include <exception>

//----------//
//   ROOT   //
//----------//   
#include <TH1.h>
#include <TMath.h>
#include <TCanvas.h>

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
#include <RooProfileLL.h>
#include <TCanvas.h>
#include <RooPlot.h>

#include <Math/ProbFunc.h>

//---//
#include <HpsFitResult.h>

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
        std::map<double, HpsFitResult*> fitWindow(TH1* histogram, double window_start, double window_end, double window_step);

        /**
         * Fit the given histogram in the window with range 
         * (window_start, window_start + window_size).  
         * 
         * @param data The RooFit histogram to fit.
         * @param window_start The start of the fit window.
         */
        HpsFitResult* fitWindow(RooDataHist* data, double window_start);

        /**
         * Fit the given histogram in the window with range 
         * (window_start, window_start + window_size).  
         * 
         * @param histogram The histogram to fit.
         * @param window_start The start of the fit window.
         */
        HpsFitResult* fitWindow(TH1* histogram, double window_start);

        
        /**
         * Fit the given histogram. If a range is specified, only fit within the 
         * range of interest.
         *
         * @param data The RooFit histogram to fit.
         * @param migrad_only If true, only run migrad.
         * @param range_name The range to fit.
         */        
        HpsFitResult* fit(RooDataHist* data, bool migrad_only, std::string range_name); 
    
        /**
         *
         */
        void calculatePValue(RooDataHist* data, HpsFitResult* result, std::string range_name, double alpha); 


        /** Fit using a background only model. */
        void fitBkgOnly();

        /** Fix the size of the window that will be fit. */
        void fixWindowSize(bool fix_window = true) { this->fix_window = fix_window; }; 

        /** Set the histogram bounds. */
        void setBounds(double low_bound, double high_bound); 

        /** Enable/disable debug */
        void setDebug(bool debug = true) { this->debug = debug; }; 

        /**
         * 
         */
        void setWindowSize(double window_size) { this->window_size = window_size; }; 

        /** Write the fit results to a text file */
        void writeResults(); 

        void getUpperLimit(TH1* histogram, HpsFitResult* result, double ap_mass);

        void getUpperLimit(RooDataHist* data, HpsFitResult* result, double ap_mass); 
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
         * Print debug statement.
         *
         * @param message Debug statement to print.
         */
        void printDebug(std::string message); 
         

        /**
         * Reset the fit parameters to their initial values.
         *
         * @param initial_params A list containing the fit parameters.
         */ 
        void resetParameters(RooArgList initial_params); 

        /**
         *
         */
        void getChi2Prob(double min_nll_null, double min_nll, double &q0, double &p_value); 

        
        /**
         *
         */
        //void generateToys(double n_toys); 

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

        /** The histogram boundary on the lower end. */
        double low_bound;

        /** The histogram boundary on the higher end. */
        double high_bound;
        
        /** Maximum size of the window */
        double max_window_size; 

        /** Size of the background window that will be used to fit. */
        double window_size;

        /** Polynomial order used to model the background. */
        int bkg_poly_order;

        /** Use a model that only includes the background. */
        bool bkg_only; 

        /** Debug flag */
        bool debug;

        /** Fix the size of the window that will be fit. */
        bool fix_window; 
};

#endif // __BUMP_HUNTER_H__
