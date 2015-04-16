/**
 * @file HpsAnalysis.h
 * @brief Interface describing an analysis
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date April 14, 2015
 *
 */

#ifndef __HPS_ANALYSIS_H__
#define __HPS_ANALYSIS_H__

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <string>

//---------------//
//--- HPS DST ---//
//---------------//
#include <HpsEvent.h>

class HpsAnalysis { 

    public: 

        /**
         *  Destructor
         */
        virtual ~Analysis() {};

        /**
         *  Method used to initialize an HPS analysis.
         */
        virtual void initialize() = 0;

        /**
         *  Method containing the code used to process an HpsEvent.
         *
         *  @param event : HpsEvent that will be processed
         */
        virtual void processEvent(HpsEvent* event) = 0;

        /**
         *  Method used to finalize an HPS analysis.
         */
        virtual void finalize() = 0;

        /**
         *  Method used to initialize any histograms used by the analysis.
         */
        // TODO:  This should use a histogram factory instead
        virtual void bookHistograms() = 0;

        /**
         *  Provide a string representation of this analysis.
         *
         *  @return String representation of this analysis.
         */
        virtual std::string toString() = 0; 
        
}; // Analysis

#endif // __ANALYSIS_H__
