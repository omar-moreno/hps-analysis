/**
 *
 * @file TridentAnalysis.h
 * @brief Analysis used to look at Tridents.
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date September 29, 2015
 *
 */

#ifndef __TRIDENT_ANALYSIS_H__
#define __TRIDENT_ANALYSIS_H__

//------------------//
//   HPS Analysis   //
//------------------//
#include <HpsAnalysis.h>
#include <Plotter.h>
#include <TrackClusterMatcher.h>
#include <EcalUtils.h>

class TridentAnalysis : public HpsAnalysis { 

    public: 

        /** Constructor */
        TridentAnalysis(); 

        /** Destructor */
        ~TridentAnalysis();

        /** Initialize an HPS analysis. */
        void initialize(); 

        /**
         * Process an HpsEvent.
         *
         * @param event HpsEvent that will be processed.
         */
        void processEvent(HpsEvent* event); 

        /** Finalize an HPS analysis. */
        void finalize(); 

        /** Initialize histograms used in this analysis. */
        void bookHistograms(); 

        /** @return A string representation of this analysis. */
        std::string toString(); 

    private: 

        /** Histogram factory used to build and save histograms. */
        Plotter* plotter;

        /** Track-Ecal cluster matcher. */
        TrackClusterMatcher* matcher; 

        /** A set of Ecal utilities */
        EcalUtils* ecal_utils;

        /** Name of the class */
        std::string class_name; 

}; // TridentAnalysis

#endif // __TRIDENT_ANALYSIS_H__
