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

//-----------------//
//--- HPS Event ---//
//-----------------//
#include <HpsEvent.h>
#include <SvtTrack.h>
#include <EcalCluster.h>

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

    private:

        /**
         * Check if an Ecal cluster and SVT track are a track match.
         *
         * @param track : The SVT track
         * @param cluster : The Ecal cluster of interest 
         */
        bool isMatch(EcalCluster* cluster, SvtTrack* track);

        /** Map from an Ecal cluster to a track matched to it */
        std::map <EcalCluster*, SvtTrack*> cluster_map;

        /** Map from an SVT track to a cluster matched to it */
        std::map <SvtTrack*, EcalCluster*> track_map;

        /** Plotter */

        double top_cluster_track_match_delta_x_low; 
        double bottom_cluster_track_match_delta_x_low; 
        double top_cluster_track_match_delta_x_high; 
        double bottom_cluster_track_match_delta_x_high; 
        
        double top_cluster_track_match_delta_y_low; 
        double bottom_cluster_track_match_delta_y_low; 
        double top_cluster_track_match_delta_y_high; 
        double bottom_cluster_track_match_delta_y_high; 
};

#endif // __TRACK_CLUSTER_MATCHER_H__
