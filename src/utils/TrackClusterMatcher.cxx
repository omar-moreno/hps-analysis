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
      // Data Moller
      top_cluster_track_match_delta_x_low(-6.10094), // data
      top_cluster_track_match_delta_x_high(12.93466), // data
      bottom_cluster_track_match_delta_x_low(-8.024475), // data
      bottom_cluster_track_match_delta_x_high(10.835475), // data
      bottom_cluster_track_match_delta_y_high(7.3967), 
      bottom_cluster_track_match_delta_y_low(-8.3093),
      top_cluster_track_match_delta_y_high(11.494),
      top_cluster_track_match_delta_y_low(-6.07562), 
      // Data FEE
      /*top_cluster_track_match_delta_x_low(-5.86232), // data
      top_cluster_track_match_delta_x_high(9.11768), // data
      bottom_cluster_track_match_delta_x_low(-5.61583), // data
      bottom_cluster_track_match_delta_x_high(8.40442), // data
      bottom_cluster_track_match_delta_y_high(6.381765), 
      bottom_cluster_track_match_delta_y_low(-7.610285),
      top_cluster_track_match_delta_y_high(7.684701),
      top_cluster_track_match_delta_y_low(-6.228299), */
      // Moller MC //
      /*top_cluster_track_match_delta_x_low(-9.4235), // moller mc
      bottom_cluster_track_match_delta_x_low(-8.6535), // moller mc
      top_cluster_track_match_delta_x_high(12.4115), // moller mc
      bottom_cluster_track_match_delta_x_high(12.2715), // moller mc
      top_cluster_track_match_delta_y_low(-8.05575), // moller mc
      bottom_cluster_track_match_delta_y_low(-13.61441), // moller mc
      top_cluster_track_match_delta_y_high(15.01425), // moller mc
      bottom_cluster_track_match_delta_y_high(7.92559), // moller mc */
      // FEE
      /*top_cluster_track_match_delta_x_low(-5.731), // fee mc
      top_cluster_track_match_delta_x_high(8.519), // fee mc
      bottom_cluster_track_match_delta_x_low(-5.86232), // fee mc
      bottom_cluster_track_match_delta_x_high(9.11768), // fee mc
      top_cluster_track_match_delta_y_low(-6.2283), // fee mc
      top_cluster_track_match_delta_y_high(7.6847), // fee mc
      bottom_cluster_track_match_delta_y_low(-7.530397), // fee mc
      bottom_cluster_track_match_delta_y_high(6.258603), // fee mc */
      enable_plots(false), 
      use_field_map(false) {
    this->bookHistograms(); 
}

TrackClusterMatcher::~TrackClusterMatcher() { 
    cluster_map.clear();
    track_map.clear();

    delete plotter;
}

void TrackClusterMatcher::findAllMatches(HpsEvent* event) { 

    //=== DEBUG
    //std::cout << "[ TrackClusterMatcher ]: Finding all track-cluster matches." << std::endl;
    //std::cout << "[ TrackClusterMatcher ]: Number of Tracks: " << event->getNumberOfTracks() << std::endl;
    //std::cout << "[ TrackClusterMatcher ]: Number of cluster : " << event->getNumberOfEcalClusters() << std::endl;
    //=== DEBUG
    
    // Clear the track and cluster maps of all previously found matches
    cluster_map.clear();
    track_map.clear();

    // Get a set of 'good' tracks from the event
    std::vector<SvtTrack*> tracks = TrackUtils::getGoodTracksList(event); 

    // Loop over all of the tracks in the event and try to find a cluster
    // match for them
    for (auto& track : tracks) { 
    
        //=== DEBUG
        //std::cout << "[ TrackClusterMatcher ]: New Track ===> Position at Ecal with field map: [ " 
        //    << track->getPositionAtEcal()[0] << ", " 
        //    << track->getPositionAtEcal()[1] << ", " 
        //    << track->getPositionAtEcal()[2] << " ]" 
        //    << std::endl;
        //=== DEBUG
        
        for (int cluster_n = 0; cluster_n < event->getNumberOfEcalClusters(); ++cluster_n) {
        
            // Get the cluster from the event 
            EcalCluster* cluster = event->getEcalCluster(cluster_n);
            
            //=== DEBUG
            //std::cout << "[ TrackClusterMatcher ]: Ecal cluster " << cluster_n << " position: [ " 
            //    << cluster->getPosition()[0] << ", " 
            //    << cluster->getPosition()[1] << ", " 
            //    << cluster->getPosition()[2] << " ]" 
            //    << std::endl;
            //=== DEBUG
            
            // Check if the track and cluster match
            double r = 0; 
            double r_min = 10000;  
            if (this->isMatch(cluster, track, r)) {
                
                //=== DEBUG
                //std::cout << "[ TrackClusterMatcher ]: Track-Cluster match found" << std::endl;
                //=== DEBUG
                
                if (r < r_min) { 
                    
                    //=== DEBUG
                    //    std::cout << "[ TrackClusterMatcher ]: New r value: " << r << std::endl;
                    //=== DEBUG
                    
                    r_min = r; 
                    cluster_map[cluster] = track;
                    track_map[track] = cluster;
                }
            }
        
            if (enable_plots && cluster_map[cluster] != nullptr) { 
                plotter->get1DHistogram("track time - cluster time - matched")->Fill(cluster->getClusterTime() 
                        - cluster_map[cluster]->getTrackTime()); 
            }
        }
    }
}

void TrackClusterMatcher::saveHistograms() { 
    
    TH2* histogram = plotter->get2DHistogram("track x @ Ecal v cluster x - track x @Ecal - top - all");
    plotter->buildGraph("average");
    TF1* gaussian = new TF1("gaussian", "gaus"); 
    for (int histogram_n = 0; histogram_n < histogram->GetNbinsX(); ++histogram_n) { 
        TH1D* projection = histogram->ProjectionY((histogram->GetName() + std::to_string(histogram_n)).c_str(), histogram_n+1, histogram_n+1);
        double mean = projection->GetMean(); 
        double rms = projection->GetRMS();
        gaussian->SetRange(-50, 50); 
        projection->Fit("gaussian", "RQ");
        plotter->add1DHistogram(projection);
        plotter->getGraph("average")->SetPoint(histogram_n, projection->GetBinCenter(histogram_n), gaussian->GetParameter(1)); 
    }
    
    delete gaussian; 

    // Save the histograms to a ROOT file
    plotter->saveToRootFile("track_cluster_matching_plots.root");
}

bool TrackClusterMatcher::isMatch(EcalCluster* cluster, SvtTrack* track, double &r) { 
  
     
    //=== DEBUG
    //std::cout << "[ TrackClusterMatcher ]: Finding match to Ecal cluster " << std::endl;
    //=== DEBUG

    // Check that the track and cluster are in the same detector volume.  If 
    // not, thre is no way they can match.
    if (track->
    if (track->isTopTrack() && cluster->getPosition()[1] < 0 
            || track->isBottomTrack() && cluster->getPosition()[1] > 0) { 
        //=== DEBUG
        //std::cout << "[ TrackClusterMatcher ]: Track and cluster are in opposite volumes" << std::endl;
        //=== DEBUG
        
        return false;
    }

    // Get the cluster position
    std::vector<double> cluster_pos = cluster->getPosition();
    
    std::vector<double> track_pos_at_ecal; 
    // Get the track position at the Ecal face.
    if (use_field_map) { 
        track_pos_at_ecal = track->getPositionAtEcal();  
    } else { 
        track_pos_at_ecal = TrackExtrapolator::extrapolateTrack(track, cluster_pos[2]);
    }

    //=== DEBUG
    //std::cout << "[ TrackClusterMatcher ]: Track position at Ecal: [ " << track_pos_at_ecal[0] << ", " 
    //    << track_pos_at_ecal[1] << ", " << track_pos_at_ecal[2] << std::endl; 
    //=== DEBUG

    double p = AnalysisUtils::getMagnitude(track->getMomentum());

    double delta_x = cluster_pos[0] - track_pos_at_ecal[0];
    double delta_y = cluster_pos[1] - track_pos_at_ecal[1];
    r = sqrt(delta_x*delta_x + delta_y*delta_y); 

    if (enable_plots) { 
        if (track->isTopTrack()) { 
           
            plotter->get1DHistogram("cluster x - track x @ Ecal - top - all")->Fill(delta_x); 
            plotter->get1DHistogram("cluster y - track y @ Ecal - top - all")->Fill(delta_y); 
            plotter->get1DHistogram("r - top - all")->Fill(r); 

            plotter->get2DHistogram("cluster x v track x @ Ecal - top - all")->Fill(cluster_pos[0], 
                    track_pos_at_ecal[0]);
            plotter->get2DHistogram("cluster x v cluster x - track x @Ecal - top - all")->Fill(cluster_pos[0], delta_x); 
            plotter->get2DHistogram("track x @ Ecal v cluster x - track x @Ecal - top - all")->Fill(track_pos_at_ecal[0],
                    delta_x); 
            plotter->get2DHistogram("cluster y v track y @ Ecal - top - all")->Fill(cluster_pos[1], 
                    track_pos_at_ecal[1]);

        } else if (track->isBottomTrack()) { 
            
            plotter->get2DHistogram("cluster x v track x @ Ecal - bottom - all")->Fill(cluster_pos[0], 
                    track_pos_at_ecal[0]);
            plotter->get2DHistogram("track x @ Ecal v cluster x - track x @Ecal - bottom - all")->Fill(track_pos_at_ecal[0],
                    delta_x); 
            plotter->get2DHistogram("cluster y v track y @ Ecal - bottom - all")->Fill(cluster_pos[1], 
                    track_pos_at_ecal[1]);
            plotter->get1DHistogram("r - bottom - all")->Fill(r); 

            plotter->get1DHistogram("cluster x - track x @ Ecal - bottom - all")->Fill(delta_x);
            plotter->get1DHistogram("cluster y - track y @ Ecal - bottom - all")->Fill(delta_y);
        }
    }

    // Check that dx and dy between the extrapolated track and cluster
    // positions is reasonable
    if ((track->isTopTrack() && (delta_x > top_cluster_track_match_delta_x_high ||
                    delta_x < top_cluster_track_match_delta_x_low)) ||
        (track->isBottomTrack() && (delta_x > bottom_cluster_track_match_delta_x_high ||
                                   delta_x < bottom_cluster_track_match_delta_x_low ))) return false;

    if (track->isTopTrack()) { 
        plotter->get1DHistogram("cluster x - track x @ Ecal - top - cuts: delta x")->Fill(delta_x);
        plotter->get1DHistogram("cluster y - track y @ Ecal - top - cuts: delta x")->Fill(delta_y);
    } else if (track->isBottomTrack()) { 
        plotter->get1DHistogram("cluster x - track x @ Ecal - bottom - cuts: delta x")->Fill(delta_x);
        plotter->get1DHistogram("cluster y - track y @ Ecal - bottom - cuts: delta x")->Fill(delta_y);
    }

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
            plotter->get2DHistogram("track x @ Ecal v cluster x - track x @Ecal - top - matched")->Fill(track_pos_at_ecal[0],
                    delta_x); 
            plotter->get2DHistogram("cluster y v track y @ Ecal - top - matched")->Fill(cluster_pos[1], 
                    track_pos_at_ecal[1]); 

            plotter->get1DHistogram("r - top - matched")->Fill(r); 

        } else if (track->isBottomTrack()) { 
            
            plotter->get2DHistogram("cluster x v track x @ Ecal - bottom - matched")->Fill(cluster_pos[0], 
                    track_pos_at_ecal[0]);
            plotter->get2DHistogram("track x @ Ecal v cluster x - track x @Ecal - bottom - matched")->Fill(track_pos_at_ecal[0],
                    delta_x); 
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
    
    TH1* plot = nullptr; 
    
    //   All tracks and clusters   //
    /////////////////////////////////
    
    //
    //   Top   
    //
    plot = plotter->build1DHistogram("cluster x - track x @ Ecal - top - all", 200, -200, 200);
    plot->GetXaxis()->SetTitle("Ecal cluster x - track x @ Ecal");

    plot = plotter->build1DHistogram("cluster y - track y @ Ecal - top - all", 100, -100, 100);
    plot->GetXaxis()->SetTitle("Ecal cluster y - track y @ Ecal");
   
    plot = plotter->build2DHistogram("track x @ Ecal : cluster x - track x @ Ecal - top - all", 200, -200, 200, 200, -200, 200);
    plot->GetYaxis()->SetTitle("Ecal cluster x - track x @ Ecal");
    
    plot = plotter->build2DHistogram("track y @ Ecal : cluster x - track x @ Ecal - top - all", 100, -100, 100, 200, -200, 200);
    plot->GetYaxis()->SetTitle("Ecal cluster x - track x @ Ecal");

    plot = plotter->build2DHistogram("track x @ Ecal : cluster y - track y @ Ecal - top - all", 200, -200, 200, 100, -100, 100);
    plot->GetYaxis()->SetTitle("Ecal cluster x - track x @ Ecal");
    
    plot = plotter->build2DHistogram("track y @ Ecal : cluster y - track y @ Ecal - top - all", 100, -100, 100, 100, -100, 100);
    plot->GetYaxis()->SetTitle("Ecal cluster x - track x @ Ecal");
   
    plot = plotter->build2DHistogram("cluster x v track x @ Ecal - top - all", 200, -200, 200, 200, -200, 200);

    plot = plotter->build2DHistogram("cluster y v track y @ Ecal - top - all", 100, -100, 100, 100, -100, 100);
    
    plot = plotter->build1DHistogram("r - top - all", 100, 0, 200);
    //
    // Bottom
    //

    plot = plotter->build1DHistogram("cluster x - track x @ Ecal - bottom - all", 200, -200, 200)
    plot->GetXaxis()->SetTitle("Ecal cluster x - track x @ Ecal");

    plot = plotter->build1DHistogram("cluster y - track y @ Ecal - bottom - all", 100, -100, 100)
    plot->GetXaxis()->SetTitle("Ecal cluster y - track y @ Ecal");

    plot = plotter->build2DHistogram("track x @ Ecal : cluster x - track x @ Ecal - bottom - all", 200, -200, 200, 200, -200, 200);
    plot->GetYaxis()->SetTitle("Ecal cluster x - track x @ Ecal");
    
    plot = plotter->build2DHistogram("track y @ Ecal : cluster x - track x @ Ecal - bottom - all", 100, -100, 100, 200, -200, 200);
    plot->GetYaxis()->SetTitle("Ecal cluster x - track x @ Ecal");

    plot = plotter->build2DHistogram("track x @ Ecal : cluster y - track y @ Ecal - bottom - all", 200, -200, 200, 100, -100, 100);
    plot->GetYaxis()->SetTitle("Ecal cluster x - track x @ Ecal");
    
    plot = plotter->build2DHistogram("track y @ Ecal : cluster y - track y @ Ecal - bottom - all", 100, -100, 100, 100, -100, 100);
    plot->GetYaxis()->SetTitle("Ecal cluster x - track x @ Ecal");


    //   After dx cuts   //
    ///////////////////////
    
    plotter->build1DHistogram("cluster x - track x @ Ecal - top - cuts: delta x", 200, -200, 200)
        ->GetXaxis()->SetTitle("Ecal cluster x - track x @ Ecal");

    plotter->build2DHistogram("cluster x v cluster x - track x @Ecal - top - all", 200, -200, 200, 200, -200, 200);
    plotter->build2DHistogram("track x @ Ecal v cluster x - track x @Ecal - top - all", 200, -200, 200, 200, -200, 200);
    //plotter->build2DHistogram("p v track x @ Ecal - top", 50, 0, 1.5, 200, -200, 200);
    //plotter->build2DHistogram("p v cluster x - top", 50, 0, 1.5, 200, -200, 200);
    //plotter->build2DHistogram("cluster pair energy v cluster x - top", 50, 0, 1.5, 200, -200, 200);
    //plotter->build2DHistogram("cluster x - track x v e/p - top", 200, -200, 200, 40, 0, 2);
   

    plotter->build1DHistogram("cluster y - track y @ Ecal - top - cuts: delta x", 100, -100, 100)
        ->GetXaxis()->SetTitle("Ecal cluster y - track y @ Ecal");

    //--- Bottom ---//
    plotter->build1DHistogram("cluster x - track x @ Ecal - bottom - cuts: delta x", 200, -200, 200)
        ->GetXaxis()->SetTitle("Ecal cluster x - track x @ Ecal");
    plotter->build2DHistogram("cluster x v track x @ Ecal - bottom - all", 200, -200, 200, 200, -200, 200);
    plotter->build2DHistogram("track x @ Ecal v cluster x - track x @Ecal - bottom - all", 200, -200, 200, 200, -200, 200);
    //plotter->build2DHistogram("p v track x @ Ecal - bottom", 50, 0, 1.5, 200, -200, 200);
    //plotter->build2DHistogram("p v cluster x - bottom", 50, 0, 1.5, 200, -200, 200);
    //plotter->build2DHistogram("cluster pair energy v cluster x - bottom", 50, 0, 1.5, 200, -200, 200);
    //plotter->build2DHistogram("cluster x - track x v e/p - bottom", 200, -200, 200, 40, 0, 2);

    plotter->build1DHistogram("cluster y - track y @ Ecal - bottom - cuts: delta x", 100, -100, 100)
        ->GetXaxis()->SetTitle("Ecal cluster y - track y @ Ecal");
    plotter->build2DHistogram("cluster y v track y @ Ecal - bottom - all", 100, -100, 100, 100, -100, 100);

    plotter->build1DHistogram("r - bottom - all", 100, 0, 200);

    //--- Matched tracks ---//
    //----------------------//
    
    //--- Top ---//
    plotter->build1DHistogram("cluster x - track x @ Ecal - top - matched", 200, -200, 200)
        ->GetXaxis()->SetTitle("Ecal cluster x - track x @ Ecal");
    plotter->build2DHistogram("cluster x v track x @ Ecal - top - matched", 200, -200, 200, 200, -200, 200);
    plotter->build2DHistogram("track x @ Ecal v cluster x - track x @Ecal - top - matched", 200, -200, 200, 200, -200, 200);
   
    plotter->build1DHistogram("cluster y - track y @ Ecal - top - matched", 100, -100, 100)
        ->GetXaxis()->SetTitle("Ecal cluster y - track y @ Ecal");
    plotter->build2DHistogram("cluster y v track y @ Ecal - top - matched", 100, -100, 100, 100, -100, 100);

    plotter->build1DHistogram("r - top - matched", 100, 0, 200);

    //--- Bottom ---//
    plotter->build1DHistogram("cluster x - track x @ Ecal - bottom - matched", 200, -200, 200)
        ->GetXaxis()->SetTitle("Ecal cluster x - track x @ Ecal");
    plotter->build2DHistogram("cluster x v track x @ Ecal - bottom - matched", 200, -200, 200, 200, -200, 200);
    plotter->build2DHistogram("track x @ Ecal v cluster x - track x @Ecal - bottom - matched", 200, -200, 200, 200, -200, 200);

    plotter->build1DHistogram("cluster y - track y @ Ecal - bottom - matched", 100, -100, 100)
        ->GetXaxis()->SetTitle("Ecal cluster y - track y @ Ecal");
    plotter->build2DHistogram("cluster y v track y @ Ecal - bottom - matched", 100, -100, 100, 100, -100, 100);

    plotter->build1DHistogram("r - bottom - matched", 100, 0, 200);

    plotter->build1DHistogram("track time - cluster time - matched", 60, 30, 60); 
}
