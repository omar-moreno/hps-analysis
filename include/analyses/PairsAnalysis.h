/**
 * @file PairsAnalysis.h
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date May 14, 2015
 *
 */

#ifndef __PAIRS_ANALYSIS_H__
#define __PAIRS_ANALYSIS_H__

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <unordered_map>

//------------//
//--- ROOT ---//
//------------//
#include <TFile.h>

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <HpsAnalysis.h>
#include <Plotter.h>

//---------------//
//--- HPS DST ---//
//---------------//
#include <HpsParticle.h>

class PairsAnalysis : public HpsAnalysis { 

    public: 

        /**
         * Constructor
         */
        PairsAnalysis();

        /**
         * Destructor
         */
        ~PairsAnalysis();

        /**
         *  Method used to initialize an HPS analysis.
         */
        void initialize();

        /**
         *  Method containing the code used to process an HpsEvent.
         *
         *  @param event : HpsEvent that will be processed
         */
        void processEvent(HpsEvent* event);

        /**
         *  Method used to finalize an HPS analysis.
         */
        void finalize();

        /**
         *  Method used to initialize any histograms used by the analysis.
         */
        // TODO:  This should use a histogram factory instead
        void bookHistograms();

        /**
         *  Provide a string representation of this analysis.
         *
         *  @return String representation of this analysis.
         */
        std::string toString();

    private:

        HpsParticle* particle;

        TCanvas* canvas; 
        
        // ROOT output_file
        TFile* output_file; 

        Plotter* uc_vtx_plotter;
        Plotter* uc_vtx_epem_plotter;
        Plotter* uc_vtx_emem_plotter;
        
        // Name of the class
        std::string class_name;

}; // PairsAnalysis

#endif // __PAIRS_ANALYSIS__
