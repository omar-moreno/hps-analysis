
#include <TrackClusterMatcher.h>

TrackClusterMatcher::TrackClusterMatcher() {
}

TrackClusterMatcher::~TrackClusterMatcher() { 
}

void TrackClusterMatcher::findAllMatches(HpsEvent* event) {

    cluster_map.clear();
    track_map.clear();    
    for (int cluster_n = 0; cluster_n < event->getNumberOfEcalClusters(); ++cluster_n) {
    
        // Get the cluster from the event 
        EcalCluster* cluster = event->getEcalCluster(cluster_n);
    
        for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
        
            // Get a track from the event
            SvtTrack* track = event->getTrack(track_n);
             
            if (this->isMatch(cluster, track)) {
                cluster_map[cluster] = track;
                track_map[track] = cluster;
                break;
            }
        } 
    } 
}


SvtTrack* TrackClusterMatcher::getMatchingTrack(EcalCluster* cluster) { 
    return cluster_map[cluster]; 
}

EcalCluster* TrackClusterMatcher::getMatchingCluster(SvtTrack* track) { 
    return track_map[track];
}

bool TrackClusterMatcher::isMatch(EcalCluster* cluster, SvtTrack* track) { 

    // Get the cluster position
    std::vector<double> cluster_pos = cluster->getPosition();

    // Extrapolate the track to the Ecal cluster position 
    std::vector<double> track_pos_at_cluster_shower_max 
        = TrackExtrapolator::extrapolateTrack(track, cluster_pos[2]);

    // If the track and cluster are in opposite volumes, then they can't 
    // be a match
    if (cluster_pos[1]*track_pos_at_cluster_shower_max[1] < 0) return false;

    // Check that dx and dy between the extrapolated track and cluster
    // positions is reasonable
    if (cluster_pos[0] - track_pos_at_cluster_shower_max[0] > 12 ||
            cluster_pos[0] - track_pos_at_cluster_shower_max[0] < -18) return false;
    
    if (cluster_pos[1] - track_pos_at_cluster_shower_max[1] > 12
            || cluster_pos[1] - track_pos_at_cluster_shower_max[1] < -18) return false;

    return true;
}


