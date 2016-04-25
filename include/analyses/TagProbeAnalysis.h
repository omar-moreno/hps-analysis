
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
#include <TrackClusterMatcher.h>
#include <Plotter.h>
#include <EcalUtils.h>

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

        bool passClusterEnergySumCut(EcalCluster* first_cluster, EcalCluster* second_cluster);

        bool passFiducialCut(EcalCluster* first_cluster, EcalCluster* second_cluster); 

        Plotter* plotter;

        /** A set of Ecal utilities */
        EcalUtils* ecal_utils;

        /** Track-Cluster matcher */
        TrackClusterMatcher* matcher;  

        /** */
        FlatTupleMaker* tuple; 
        

        //   Event counters   //
        //--------------------//

        /** Total number of triggers */
        double event_count;
        
        /** */ 
        double good_ecal_pair_count; 

        /** */
        double fiducial_cut_pass_count;


        double ecal_cluster_sum_cut_pass_count;
        double ecal_cluster_diff_cut_pass_count;
        double ecal_cluster_time_cut_pass_count;
        double ecal_row1_cut_pass_count;
        double ecal_cluster_y_cut_pass_count;  
        double tag_ep_cut_pass_count; 
        double top_tag_cluster_count; 
        double bottom_tag_cluster_count; 
        double top_tag_cand_cluster_count; 
        double bottom_tag_cand_cluster_count; 

        // Name of the class
		std::string class_name;


};  // TagProbeAnalysis

#endif // __TAG_PROBE_ANALYSIS_H__
