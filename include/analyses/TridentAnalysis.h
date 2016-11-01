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
#include <FlatTupleMaker.h>
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
    
    protected: 

        /** Name of the class */
        std::string class_name; 
    
    private: 

        /** Utility used to create ROOT ntuples. */
        FlatTupleMaker* tuple{new FlatTupleMaker("trident_analysis.root", "results")}; 

        /** Track-Ecal cluster matcher. */
        TrackClusterMatcher* matcher; 

        /** A set of Ecal utilities */
        EcalUtils* ecal_utils;

        /** Total number of events processed */
        double event_counter{0};

        /** Total number of events with tracks */
        double event_has_track{0};

        /** Total number of events with a positron */
        double event_has_positron{0};

        double good_cluster_pair_counter; 
        double matched_event_counter; 
        double v0_cand_counter;

}; // TridentAnalysis

#endif // __TRIDENT_ANALYSIS_H__
