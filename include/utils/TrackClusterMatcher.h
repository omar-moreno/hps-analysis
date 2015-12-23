/**
 * @file TrackClusterMatcher.h
 * @brief  A class used to find all Ecal cluster - track matches
 * @author Omar Moreno <omoreno1@ucsc.edu> \n
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date June 16, 2015 
 */


#ifndef __TRACK_CLUSTER_MATCHER_H__
#define __TRACK_CLUSTER_MATCHER_H__

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <map>

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <TrackExtrapolator.h>
#include <Plotter.h>
#include <TrackType.h>
#include <StrategyType.h>
#include <TrackUtils.h>

//-----------------//
//--- HPS Event ---//
//-----------------//
#include <HpsEvent.h>
#include <SvtTrack.h>
#include <EcalCluster.h>

#include <TF1.h>

class TrackClusterMatcher {

    public: 

        /**
         * Constructor.
         */
        TrackClusterMatcher();

        /**
         * Destructor.
         */
        ~TrackClusterMatcher();

        /**
         * Find all cluster and track matches in the event.
         *
         * @param event : HpsEvent containing the tracks and clusters of 
         *                interest
         */
        void findAllMatches(HpsEvent* event);

        /**
         * Get a track that is matched to the given cluster.
         *
         * @param cluster : The Ecal cluster of interest 
         */
        SvtTrack* getMatchingTrack(EcalCluster* cluster) { return cluster_map[cluster]; };

        /**
         * Get a cluster that is matched to the given track.
         *
         * @param track : The SVT track
         */
        EcalCluster* getMatchingCluster(SvtTrack* track) { return track_map[track]; };

        /**
         * Use the extrapolated track position at the Ecal face found using
         * the full field map.
         *
         * @param use_field_map If true, use the extrapolated track position 
         *                      found using the field map.  If false, 
         *                      extrapolated the track analytically.
         */
        void useFieldMap(bool use_field_map) { this->use_field_map = use_field_map; }; 

        /**
         * Enable/disable booking, filling and saving of plots.
         *
         * @param enable_plots : true to enable, false to disable
         */
        void enablePlots(bool enable_plots = true) { this->enable_plots = enable_plots; }; 

        /**
         * Save the histograms to a ROOT file.
         */
        void saveHistograms(); 

    private:

        /**
         * Check if an Ecal cluster and SVT track are a track match.
         *
         * @param track : The SVT track
         * @param cluster : The Ecal cluster of interest 
         */
        bool isMatch(EcalCluster* cluster, SvtTrack* track, double &r);

        /**
         * Book histograms.
         */
        void bookHistograms(); 

        /** Map from an Ecal cluster to a track matched to it */
        std::map <EcalCluster*, SvtTrack*> cluster_map;

        /** Map from an SVT track to a cluster matched to it */
        std::map <SvtTrack*, EcalCluster*> track_map;

        /** Plotter */
        Plotter* plotter; 

        double top_cluster_track_match_delta_x_low; 
        double bottom_cluster_track_match_delta_x_low; 
        double top_cluster_track_match_delta_x_high; 
        double bottom_cluster_track_match_delta_x_high; 
        
        double top_cluster_track_match_delta_y_low; 
        double bottom_cluster_track_match_delta_y_low; 
        double top_cluster_track_match_delta_y_high; 
        double bottom_cluster_track_match_delta_y_high;

        bool enable_plots; 
        
        /** 
         * Flag indicating whether the extrapolated track position at the Ecal
         * face found using the full field map should be used. 
         */
        bool use_field_map; 
};

#endif // __TRACK_CLUSTER_MATCHER_H__
