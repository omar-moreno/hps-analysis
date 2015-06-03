
#include <TagProbeAnalysis.h>

TagProbeAnalysis::TagProbeAnalysis() 
    : track(NULL),
      matched_tag_track(NULL),
      cluster(NULL),
      ecal_plotter(NULL),
      svt_plotter(NULL),
      cluster_energy_threshold(.200),
      cluster_energy_sum_threshold(.800),
      total_events(0),
      pass_fiducial_cut(0),  
      pass_cluster_threshold_cut(0),
      pass_cluster_energy_sum_cut(0),  
      candidates(0),
      found(0), 
      class_name("TagProbeAnalysis") {       
}

TagProbeAnalysis::~TagProbeAnalysis() { 
}

void TagProbeAnalysis::initialize() { 
    this->bookHistograms();
}

void TagProbeAnalysis::processEvent(HpsEvent* event) {

    if (!event->isPair1Trigger()) return;

    // Only look at events that have two Ecal clusters
    if (event->getNumberOfEcalClusters() != 2) return;
    ++total_events;

    if (!passFiducialCut(event)) return;
    ++pass_fiducial_cut;

    // Check that both clusters pass the cluster energy cut
    if (!passClusterEnergyCut(event)) return;
    ++pass_cluster_threshold_cut;

    // Check that the sum of the cluster energy is less than the beam energy
    if (!passClusterEnergySumCut(event)) return;
    ++pass_cluster_energy_sum_cut;  

    // Randomly choose one of the Ecal clusters.
    srand(time(0)); 
    int cluster_index = rand()%2;
    cluster = event->getEcalCluster(cluster_index);

    //std::cout << "[ TagProbeAnalysis ]: Number of tracks: " << event->getNumberOfTracks() << std::endl;
    matched_tag_track = NULL;
    
    // Try and match a track to the cluster.
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
        
        // Get a track from the event
        track = event->getTrack(track_n);
    
        if (isMatch(cluster, track)) { 
            matched_tag_track = track;
            break;
        }
    }

    if (matched_tag_track == NULL) return; 

    candidates++;

    if (cluster_index == 0) { 
        cluster = event->getEcalCluster(1);
    } else {
        cluster = event->getEcalCluster(0);
    }

    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
        
        // Get a track from the event
        track = event->getTrack(track_n);
       
        if (matched_tag_track == track) continue;  

        if (isMatch(cluster, track)) { 
            found++;
        }
    } 
}

void TagProbeAnalysis::finalize() { 
    
    std::cout << "//---------------------------------------------------//" << std::endl;
    std::cout << "// Total events: " << total_events << std::endl;
    std::cout << "// Events passing cluster threshold cut: " 
        << pass_cluster_threshold_cut/total_events << std::endl;
    std::cout << "// Events passing cluster energy sum cut: " 
        << pass_cluster_energy_sum_cut/total_events << std::endl;

    std::cout << "Candidates: " << candidates << std::endl;
    std::cout << "Found: " << found << std::endl;

    std::cout << "Tracking efficiency: " << (found/candidates)*100 << std::endl;
    std::cout << "//---------------------------------------------------//" << std::endl;

    return;
}

void TagProbeAnalysis::bookHistograms() { 
    
    ecal_plotter->setType("float");
    ecal_plotter->build1DHistogram("Cluster energy - All", 50, 0, 2);
    ecal_plotter->build1DHistogram("Cluster energy - Pass Cut", 50, 0, 2);
    ecal_plotter->build2DHistogram("Probe clusters", 47, -23, 24, 12, -6, 6);

    svt_plotter->setType("float");
    svt_plotter->build1DHistogram("Track-Cluster dx", 100, -200, 200);
    svt_plotter->build1DHistogram("Track-Cluster dy", 100, -200, 200);

    return;
}

std::string TagProbeAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;
}

bool TagProbeAnalysis::passClusterEnergyCut(HpsEvent* event) { 

    for (int cluster_n = 0; cluster_n < event->getNumberOfEcalClusters(); ++cluster_n) { 
        if (event->getEcalCluster(cluster_n)->getEnergy() < cluster_energy_threshold) return false; 
    }
    return true;
}

bool TagProbeAnalysis::passClusterEnergySumCut(HpsEvent* event) {
    
    double cluster_energy_sum = 0; 
    for (int cluster_n = 0; cluster_n < event->getNumberOfEcalClusters(); ++cluster_n) { 
        cluster_energy_sum += event->getEcalCluster(cluster_n)->getEnergy();
    }

    if (cluster_energy_sum > cluster_energy_sum_threshold) return false; 

    return true;
}

bool TagProbeAnalysis::passFiducialCut(HpsEvent* event) { 
    
    if (event->getEcalCluster(0)->getPosition()[1]*event->getEcalCluster(1)->getPosition()[1] > 0) return false;

    return true;
}

bool TagProbeAnalysis::isMatch(EcalCluster* cluster, SvtTrack* track) { 

    // Get the cluster position
    std::vector<double> cluster_pos = cluster->getPosition();
    /*std::cout << "[ TagProbeAnalysis ]: ECal cluster position: " 
        << " x: " << cluster_pos[0] 
        << " y: " << cluster_pos[1] 
        << " z: " << cluster_pos[2]
        << std::endl;*/

    // Extrapolate the track to the Ecal cluster position 
    std::vector<double> track_pos_at_cluster_shower_max 
        = TrackExtrapolator::extrapolateTrack(track, cluster_pos[2]);
    /*std::cout << "[ TagProbeAnalysis ]: Track position at shower max: " 
         << " x: " << track_pos_at_cluster_shower_max[0] 
         << " y: " << track_pos_at_cluster_shower_max[1] 
         << " z: " << track_pos_at_cluster_shower_max[2]
         << std::endl;*/ 

    // If the track and cluster are in opposite volumes, then they can't 
    // be a match
    if (cluster_pos[1]*track_pos_at_cluster_shower_max[1] < 0) return false;
   
    // Check that dx and dy between the extrapolated track and cluster
    // positions is reasonable
    if (abs(cluster_pos[0] - track_pos_at_cluster_shower_max[0]) > 50) return false;
    
    if (abs(cluster_pos[1] - track_pos_at_cluster_shower_max[1]) > 50) return false;

    return true;
}
