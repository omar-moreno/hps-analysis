
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

    // Pick one of the two clusters at random then try and match
    // a track to it.
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
        
        // Get a track from the event
        track = event->getTrack(track_n);

        // Check that the track and the cluster are in the same detector volume
        if (!(track->isTopTrack() && cluster_position[1] > 0) || 
                !(track->isBottomTrack() && cluster_position[1] < 0)) continue;

        // Extrapolate the track to the Ecal cluster position 
        
         
    }
}
