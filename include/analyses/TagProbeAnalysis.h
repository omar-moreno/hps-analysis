
#ifndef __TAG_PROBE_ANALYSIS_H__
#define __TAG_PROBE_ANALYSIS_H__

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <vector>
#include <ctime>
#include <cstdlib>

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <HpsAnalysis.h>
#include <TrackExtrapolator.h>
#include <Plotter.h>

//---------------//
//--- HPS DST ---//
//---------------//
#include <SvtTrack.h>
#include <EcalCluster.h>


class TagProbeAnalysis : public HpsAnalysis { 

    public: 

        /**
         * Constructor
         */
        TagProbeAnalysis(); 

        /**
         * Destructor
         */
        ~TagProbeAnalysis();

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

        bool passClusterEnergyCut(HpsEvent* event);

        bool passClusterEnergySumCut(HpsEvent* event);

        bool passFiducialCut(HpsEvent* event); 

        bool passClusterTimeCut(HpsEvent* event);  

        bool isMatch(EcalCluster* cluster, SvtTrack* track);

        Plotter* plotter;

        double total_events;
        double total_pair_trigger_events;
        double total_pair_events;
        double total_tag_candidates;

        // Name of the class
		std::string class_name;


};  // TagProbeAnalysis

#endif // __TAG_PROBE_ANALYSIS_H__
