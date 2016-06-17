
#ifndef __FEE_ANALYSIS_H__
#define __FEE_ANALYSIS_H__

#include <algorithm>

#include <TGraphAsymmErrors.h>

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <HpsAnalysis.h>
#include <TrackExtrapolator.h>
#include <Plotter.h>
#include <TrackClusterMatcher.h>
#include <RooFitter.h>

//---------------//
//--- HPS DST ---//
//---------------//
#include <SvtTrack.h>
#include <EcalCluster.h>
#include <GblTrack.h>

#include <RooRealVar.h>
#include <RooPlot.h>

#include <TFile.h>

class FeeAnalysis : public HpsAnalysis {

    public:

        /** Constructor */
        FeeAnalysis();

        /** Destructor */
        ~FeeAnalysis();

        /** 
         * Method used to initialize an HPS analysis. This method is called
         * once at the beginning of an analysis. 
         */
        void initialize();

        /**
         *  Method used to process an HpsEvent.
         *
         *  @param event : HpsEvent that will be processed
         */
        void processEvent(HpsEvent* event);

        /** 
         * Method used to finalize an HPS analysis. This method is called once
         * at the end of an analysis.
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
        bool isEdgeCrystal(EcalHit* hit);

        EcalCluster* cluster;

        TrackClusterMatcher* matcher; 

        Plotter* plotter;
        
        double cluster_energy_low_threshold; 
        double cluster_energy_high_threshold; 

        bool cuts_enabled; 

        // Name of the class
		std::string class_name;
       
        //-- Event counters --//
        //--------------------//
        
        /** Total number of trigger events. */
        int trigger_count;
        
        /** Total number of singles1 trigger events that had tracks */
        int event_track_counter;
        
        /** Total number of singles1 trigger events where the SVT bias was on */
        int bias_on_counter; 

        /** Total number of singles1 triggers */
        int single1_trigger_counter;

        int ecal_cluster_time_cut_pass_counter; 
        int ecal_cluster_size_cut_pass_counter; 
        int ecal_cluster_energy_cut_pass_counter; 
        int ecal_cluster_seed_energy_cut_pass_counter; 

        /** 
         * Total number of singles 1 triggers where the SVT bias was on and it
         * was in closed position.
         */
        int svt_closed_position_counter; 

};

#endif // __SIMPLE_TRACKING_EFFICIENCY_ANALYSIS_H__
