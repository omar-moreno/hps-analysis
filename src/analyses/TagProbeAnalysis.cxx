
#include <TagProbeAnalysis.h>

TagProbeAnalysis::TagProbeAnalysis() 
    : track(NULL),
      cluster(NULL),
      class_name("TagProbeAnalysis") {       
}

TagProbeAnalysis::~TagProbeAnalysis() { 
}

void TagProbeAnalysis::initialize() { 
    this->bookHistograms();
}

void TagProbeAnalysis::processEvent(HpsEvent* event) {

    // Only look at events that have two Ecal clusters
    if (event->getNumberOfEcalClusters() != 2) return;

    cluster = event->getEcalCluster(0);

    //
    std::vector<double> cluster_position = cluster->getPosition();
    std::cout << "[ TagProbeAnalysis ]: ECal cluster position: " 
        << " x: " << cluster_position[0] 
        << " y: " << cluster_position[1] 
        << " z: " << cluster_position[2]
        << std::endl; 

    std::cout << "[ TagProbeAnalysis ]: Number of tracks: " << event->getNumberOfTracks() << std::endl;
    double r_max = 10000;  
    // Pick one of the two clusters at random then try and match
    // a track to it.
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
        
        // Get a track from the event
        track = event->getTrack(track_n);

        // Check that the track and the cluster are in the same detector volume
        //if (!(track->isTopTrack() && cluster_position[1] > 0) || 
        //        !(track->isBottomTrack() && cluster_position[1] < 0)) continue;

        // Extrapolate the track to the Ecal cluster position 
        std::vector<double> track_position_at_shower_max = TrackExtrapolator::extrapolateTrack(track, cluster_position[2]);
        
        std::cout << "[ TagProbeAnalysis ]: Track position at shower max: " 
            << " x: " << track_position_at_shower_max[0] 
            << " y: " << track_position_at_shower_max[1] 
            << " z: " << track_position_at_shower_max[2]
            << std::endl; 
         
        double r_mag = 0;
        for (int index = 0; index < 3; ++index) { 
            double r_i = track_position_at_shower_max[index] - cluster_position[index];
            r_mag += r_i*r_i;
        }
        r_mag = sqrt(r_mag);

        std::cout << "[ TagProbeAnalysis ]: Distance between track and cluster: " << r_mag << std::endl;

        if (r_mag < r_max) { 
            r_max = r_mag; 
        } 
    }
}

void TagProbeAnalysis::finalize() { 
    return;
}

void TagProbeAnalysis::bookHistograms() { 
    return;
}

std::string TagProbeAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;
}
