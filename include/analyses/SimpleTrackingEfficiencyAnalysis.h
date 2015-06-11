
#ifndef __SIMPLE_TRACKING_EFFICIENCY_ANALYSIS_H__
#define __SIMPLE_TRACKING_EFFICIENCY_ANALYSIS_H__

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

class SimpleTrackingEfficiencyAnalysis : public HpsAnalysis {

    public:

        /**
         * Constructor
         */
        SimpleTrackingEfficiencyAnalysis();

        /**
         * Destructor
         */
        ~SimpleTrackingEfficiencyAnalysis();

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

        //--- Cuts ---//
        bool passEnergyCut(EcalCluster* cluster);
        bool passClusterTimeCut(EcalCluster* cluster);
        bool passClusterSizeCut(EcalCluster* cluster);

        bool isMatch(EcalCluster* cluster, SvtTrack* track);

        SvtTrack* track;
        EcalCluster* cluster;

        Plotter* plotter;
        
        double cluster_energy_low_threshold; 
        double cluster_energy_high_threshold; 

        bool cuts_enabled; 

        // Name of the class
		std::string class_name;
};

#endif // __SIMPLE_TRACKING_EFFICIENCY_ANALYSIS_H__
