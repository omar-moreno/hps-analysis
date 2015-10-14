/**
 *
 * @file HpsAnalysis.h
 * @brief Interface outlining an HPS analysis.
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date April 14, 2015
 *
 */

#ifndef __HPS_ANALYSIS_H__
#define __HPS_ANALYSIS_H__

//----------------//
//   C++ StdLib   //
//----------------//
#include <string>

//-------------//
//   HPS DST   //
//-------------//
#include <HpsEvent.h>

class HpsAnalysis { 

    public: 

        /** Destructor */
        virtual ~HpsAnalysis() {};

        /** 
         * Initialize an HPS analysis.  This method is called only once before
         * any events have been processed.
         */
        virtual void initialize() = 0;

        /**
         * Process an HPS event i.e. {@link HpsEvent} object.  All analysis
         * code will usually be placed within this method.
         *
         * @param event {@link HpsEvent} object to process.
         */
        virtual void processEvent(HpsEvent* event) = 0;

        /**  
         * Finalize an HPS analysis.  This method is called only once after all
         * events have been processed.
         */
        virtual void finalize() = 0;

        /**
         * Instantiate (book) any histograms that will be used by the analysis.
         * This method is here for organizational reasons as most of this code
         * could be placed within the initialize method.
         */
        virtual void bookHistograms() = 0;

        /**
         * Provide a string representation of this analysis.
         *
         * @return String representation of this analysis.
         */
        virtual std::string toString() = 0; 
        
}; // Analysis

#endif // __ANALYSIS_H__
