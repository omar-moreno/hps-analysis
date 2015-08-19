/**
 * @file TrackClusterMatcher.h
 * @brief  A class used to find all Ecal cluster - track matches
 * @author Omar Moreno <omoreno1@ucsc.edu> \n
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date June 16, 2015 
 */

#include <TrackClusterMatcher.h>

TrackClusterMatcher::TrackClusterMatcher() 
    : top_cluster_track_match_delta_x_low(-7.61),
      bottom_cluster_track_match_delta_x_low(-13.75),
      top_cluster_track_match_delta_x_high(12),
      bottom_cluster_track_match_delta_x_high(6),
      top_cluster_track_match_delta_y_low(-14),
      bottom_cluster_track_match_delta_y_low(-14),
      top_cluster_track_match_delta_y_high(14),
      bottom_cluster_track_match_delta_y_high(14)
{
}

TrackClusterMatcher::~TrackClusterMatcher() { 
    cluster_map.clear();
    track_map.clear(); 
}

void TrackClusterMatcher::findAllMatches(HpsEvent* event) {

    // Clear the track and cluster maps of all previously found track and 
    // cluster matches.
    cluster_map.clear();
    track_map.clear();

    // Loop over all of the clusters in the event and try to find a track
    // match for it.
    for (int cluster_n = 0; cluster_n < event->getNumberOfEcalClusters(); ++cluster_n) {
    
        // Get the cluster from the event 
        EcalCluster* cluster = event->getEcalCluster(cluster_n);
    
        // Loop over all of the tracks in the event and try to find a match to
        // a cluster
        for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
        
            // Get a track from the event
            SvtTrack* track = event->getTrack(track_n);
            
            // Check if the track and cluster match 
            if (this->isMatch(cluster, track)) {
                cluster_map[cluster] = track;
                track_map[track] = cluster;
                break;
            }
        } 
    } 
}

bool TrackClusterMatcher::isMatch(EcalCluster* cluster, SvtTrack* track) { 
    
    // Check that the track and cluster are in the same detector volume.
    // If not, thre is no way they can match.
    if (track->isTopTrack() && cluster->getPosition()[1] < 0
            || track->isBottomTrack() && cluster->getPosition()[1] > 0) return false;

    // Get the cluster position
    std::vector<double> cluster_pos = cluster->getPosition();
    
    // Extrapolate the track to the Ecal cluster position 
    std::vector<double> track_pos_at_ecal 
        = TrackExtrapolator::extrapolateTrack(track, cluster_pos[2]);
    
    double p = AnalysisUtils::getMagnitude(track->getMomentum());

    double delta_x = cluster_pos[0] - track_pos_at_ecal[0];
    double delta_y = cluster_pos[1] - track_pos_at_ecal[1];

    // Check that dx and dy between the extrapolated track and cluster
    // positions is reasonable
    if ((track->isTopTrack() && (delta_x > top_cluster_track_match_delta_x_high ||
                    delta_x < top_cluster_track_match_delta_x_low)) ||
        (track->isBottomTrack() && (delta_x > bottom_cluster_track_match_delta_x_high ||
                                   delta_x < bottom_cluster_track_match_delta_x_low ))) return false;
    
    if ((track->isTopTrack() && (delta_y > top_cluster_track_match_delta_y_high ||
                    delta_y < top_cluster_track_match_delta_y_low)) ||
        (track->isBottomTrack() && (delta_y > bottom_cluster_track_match_delta_y_high ||
                                    delta_y < bottom_cluster_track_match_delta_y_low))) return false;

    return true;
}


