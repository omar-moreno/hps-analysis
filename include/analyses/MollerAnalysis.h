/**
 * @file MollerAnalysis.h
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date July 10, 2015
 */

#ifndef __MOLLER_ANALYSIS_H__
#define __MOLLER_ANALYSIS_H__

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <HpsAnalysis.h>
#include <Plotter.h>
#include <AnalysisUtils.h>
#include <TrackExtrapolator.h>

//---------------//
//--- HPS DST ---//
//---------------//
#include <SvtTrack.h>

class MollerAnalysis : public HpsAnalysis { 

    public: 

        /**
         * Constructor
         */
        MollerAnalysis();

        /**
         * Destructor
         */
        ~MollerAnalysis();

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

        Plotter* plotter;
        
        // Name of the class
        std::string class_name;

};

#endif // __MOLLER_ANALYSIS_H__
