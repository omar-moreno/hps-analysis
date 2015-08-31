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
    : plotter(new Plotter()), 
      top_cluster_track_match_delta_x_low(-7.61),
      bottom_cluster_track_match_delta_x_low(-13.75),
      top_cluster_track_match_delta_x_high(12),
      bottom_cluster_track_match_delta_x_high(6),
      top_cluster_track_match_delta_y_low(-14),
      bottom_cluster_track_match_delta_y_low(-14),
      top_cluster_track_match_delta_y_high(14),
      bottom_cluster_track_match_delta_y_high(14), 
      enable_plots(false)
{
    this->bookHistograms(); 
}

TrackClusterMatcher::~TrackClusterMatcher() { 
    cluster_map.clear();
    track_map.clear();

    delete plotter;
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

void TrackClusterMatcher::saveHistograms() { 
    
    // Save the histograms to a ROOT file
    plotter->saveToRootFile("track_cluster_matching_plots.root");
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
    double r = sqrt(delta_x*delta_x + delta_y*delta_y); 

    if (enable_plots) { 
        if (track->isTopTrack()) { 
           
            plotter->get1DHistogram("cluster x - track x @ Ecal - top - all")->Fill(delta_x); 
            plotter->get1DHistogram("cluster y - track y @ Ecal - top - all")->Fill(delta_y); 
            plotter->get1DHistogram("r - top - all")->Fill(r); 

            plotter->get2DHistogram("cluster x v track x @ Ecal - top - all")->Fill(cluster_pos[0], 
                    track_pos_at_ecal[0]);
            plotter->get2DHistogram("cluster y v track y @ Ecal - top - all")->Fill(cluster_pos[1], 
                    track_pos_at_ecal[1]); 

        } else if (track->isBottomTrack()) { 
        
            plotter->get1DHistogram("cluster x - track x @ Ecal - bottom - all")->Fill(delta_x);
            plotter->get1DHistogram("cluster y - track y @ Ecal - bottom - all")->Fill(delta_y);
            plotter->get1DHistogram("r - bottom - all")->Fill(r); 

            plotter->get2DHistogram("cluster x v track x @ Ecal - bottom - all")->Fill(cluster_pos[0], 
                    track_pos_at_ecal[0]);
            plotter->get2DHistogram("cluster y v track y @ Ecal - bottom - all")->Fill(cluster_pos[1], 
                    track_pos_at_ecal[1]);
        }
    }

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

    if (enable_plots) { 
        if (track->isTopTrack()) { 
           
            plotter->get1DHistogram("cluster x - track x @ Ecal - top - matched")->Fill(delta_x); 
            plotter->get1DHistogram("cluster y - track y @ Ecal - top - matched")->Fill(delta_y); 

            plotter->get2DHistogram("cluster x v track x @ Ecal - top - matched")->Fill(cluster_pos[0], 
                    track_pos_at_ecal[0]);
            plotter->get2DHistogram("cluster y v track y @ Ecal - top - matched")->Fill(cluster_pos[1], 
                    track_pos_at_ecal[1]); 

            plotter->get1DHistogram("r - top - matched")->Fill(r); 

        } else if (track->isBottomTrack()) { 
            
            plotter->get2DHistogram("cluster x v track x @ Ecal - bottom - matched")->Fill(cluster_pos[0], 
                    track_pos_at_ecal[0]);
            plotter->get2DHistogram("cluster y v track y @ Ecal - bottom - matched")->Fill(cluster_pos[1], 
                    track_pos_at_ecal[1]);

            plotter->get1DHistogram("cluster x - track x @ Ecal - bottom - matched")->Fill(delta_x);
            plotter->get1DHistogram("cluster y - track y @ Ecal - bottom - matched")->Fill(delta_y);

            plotter->get1DHistogram("r - bottom - matched")->Fill(r); 
        }
    }

    return true;
}

void TrackClusterMatcher::bookHistograms() { 
    
    //--- All tracks and clusters ---//
    //-------------------------------//
    
    //--- Top ---//
    plotter->build1DHistogram("cluster x - track x @ Ecal - top - all", 200, -200, 200)
        ->GetXaxis()->SetTitle("Ecal cluster x - track x @ Ecal");
    plotter->build2DHistogram("cluster x v track x @ Ecal - top - all", 200, -200, 200, 200, -200, 200);
    //plotter->build2DHistogram("p v track x @ Ecal - top", 50, 0, 1.5, 200, -200, 200);
    //plotter->build2DHistogram("p v cluster x - top", 50, 0, 1.5, 200, -200, 200);
    //plotter->build2DHistogram("cluster pair energy v cluster x - top", 50, 0, 1.5, 200, -200, 200);
    //plotter->build2DHistogram("cluster x - track x v e/p - top", 200, -200, 200, 40, 0, 2);
   
    plotter->build1DHistogram("cluster y - track y @ Ecal - top - all", 100, -100, 100)
        ->GetXaxis()->SetTitle("Ecal cluster y - track y @ Ecal");
    plotter->build2DHistogram("cluster y v track y @ Ecal - top - all", 100, -100, 100, 100, -100, 100);

    plotter->build1DHistogram("r - top - all", 100, 0, 200);

    //--- Bottom ---//
    plotter->build1DHistogram("cluster x - track x @ Ecal - bottom - all", 200, -200, 200)
        ->GetXaxis()->SetTitle("Ecal cluster x - track x @ Ecal");
    plotter->build2DHistogram("cluster x v track x @ Ecal - bottom - all", 200, -200, 200, 200, -200, 200);
    //plotter->build2DHistogram("p v track x @ Ecal - bottom", 50, 0, 1.5, 200, -200, 200);
    //plotter->build2DHistogram("p v cluster x - bottom", 50, 0, 1.5, 200, -200, 200);
    //plotter->build2DHistogram("cluster pair energy v cluster x - bottom", 50, 0, 1.5, 200, -200, 200);
    //plotter->build2DHistogram("cluster x - track x v e/p - bottom", 200, -200, 200, 40, 0, 2);

    plotter->build1DHistogram("cluster y - track y @ Ecal - bottom - all", 100, -100, 100)
        ->GetXaxis()->SetTitle("Ecal cluster y - track y @ Ecal");
    plotter->build2DHistogram("cluster y v track y @ Ecal - bottom - all", 100, -100, 100, 100, -100, 100);

    plotter->build1DHistogram("r - bottom - all", 100, 0, 200);

    //--- Matched tracks ---//
    //----------------------//
    
    //--- Top ---//
    plotter->build1DHistogram("cluster x - track x @ Ecal - top - matched", 200, -200, 200)
        ->GetXaxis()->SetTitle("Ecal cluster x - track x @ Ecal");
    plotter->build2DHistogram("cluster x v track x @ Ecal - top - matched", 200, -200, 200, 200, -200, 200);
   
    plotter->build1DHistogram("cluster y - track y @ Ecal - top - matched", 100, -100, 100)
        ->GetXaxis()->SetTitle("Ecal cluster y - track y @ Ecal");
    plotter->build2DHistogram("cluster y v track y @ Ecal - top - matched", 100, -100, 100, 100, -100, 100);

    plotter->build1DHistogram("r - top - matched", 100, 0, 200);

    //--- Bottom ---//
    plotter->build1DHistogram("cluster x - track x @ Ecal - bottom - matched", 200, -200, 200)
        ->GetXaxis()->SetTitle("Ecal cluster x - track x @ Ecal");
    plotter->build2DHistogram("cluster x v track x @ Ecal - bottom - matched", 200, -200, 200, 200, -200, 200);

    plotter->build1DHistogram("cluster y - track y @ Ecal - bottom - matched", 100, -100, 100)
        ->GetXaxis()->SetTitle("Ecal cluster y - track y @ Ecal");
    plotter->build2DHistogram("cluster y v track y @ Ecal - bottom - matched", 100, -100, 100, 100, -100, 100);

    plotter->build1DHistogram("r - bottom - matched", 100, 0, 200);
}
