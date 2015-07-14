
#ifndef __TRACK_CLUSTER_MATCHER_H__
#define __TRACK_CLUSTER_MATCHER_H__

#include <map>

#include <TrackExtrapolator.h>

#include <HpsEvent.h>
#include <SvtTrack.h>
#include <EcalCluster.h>

class TrackClusterMatcher {

    public: 

        /**
         * Constructor
         */
        TrackClusterMatcher();

        /**
         *
         */
        ~TrackClusterMatcher();

        /**
         *
         */
        void findAllMatches(HpsEvent* event);

        /**
         *
         */
        SvtTrack* getMatchingTrack(EcalCluster* cluster);

        /**
         *
         */
        EcalCluster* getMatchingCluster(SvtTrack* track);

    private:

        bool isMatch(EcalCluster* cluster, SvtTrack* track);

        std::map <EcalCluster*, SvtTrack*> cluster_map;
        std::map <SvtTrack*, EcalCluster*> track_map;
};

#endif // __TRACK_CLUSTER_MATCHER_H__
