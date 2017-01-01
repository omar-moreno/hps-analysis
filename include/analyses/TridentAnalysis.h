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

        std::map<GblTrack*, int> buildSharedHitMap(HpsEvent* event);

        bool electronsShareHits(std::vector<HpsParticle*> particles, std::map<GblTrack*, int> shared_hit_map); 

        HpsParticle* getBestElectronChi2(std::vector<HpsParticle*> particles);

        bool passFeeCut(HpsParticle* particle); 

        void printDebug(std::string message); 

        /** Utility used to create ROOT ntuples. */
        FlatTupleMaker* tuple{new FlatTupleMaker("trident_analysis.root", "results")}; 

        /** Track-Ecal cluster matcher. */
        TrackClusterMatcher* matcher{new TrackClusterMatcher()}; 

        /** A set of Ecal utilities */
        EcalUtils* ecal_utils{new EcalUtils()};

        /** Total number of events processed */
        double event_counter{0};

        /** Total number of events with tracks */
        double event_has_track{0};

        /** Total number of events with a positron */
        double event_has_positron{0};

        /** Total number of events with only a single positron. */
        double event_has_single_positron{0};

        /** Total number of events with a good cluster pair */
        double event_has_good_cluster_pair{0};

        double _dumped_positron_count{0}; 

        double total_v0_good_cluster_pair{0}; 

        double total_v0_good_track_match{0};

        double _total_v0_pass_fee{0};

        bool _debug{false};

}; // TridentAnalysis

#endif // __TRIDENT_ANALYSIS_H__
