
#include <SimpleTrackingEfficiencyAnalysis.h>

SimpleTrackingEfficiencyAnalysis::SimpleTrackingEfficiencyAnalysis() 
    : track(NULL),
      cluster(NULL),
      plotter(new Plotter()),
      cluster_energy_low_threshold(.8 /* GeV */),
      cluster_energy_high_threshold(1.15 /* GeV */),
      cuts_enabled(true),
      class_name("SimpleTrackingEfficiencyAnalysis") {

}

SimpleTrackingEfficiencyAnalysis::~SimpleTrackingEfficiencyAnalysis() {
}

void SimpleTrackingEfficiencyAnalysis::initialize() { 
    this->bookHistograms();
}

void SimpleTrackingEfficiencyAnalysis::processEvent(HpsEvent* event) { 
  
    // Only look at single 0 triggers
    if (!event->isSingle0Trigger()) return;

    // For now, only look at events with a single cluster    
    if (event->getNumberOfEcalClusters() != 1) return;

    // Loop over all clusters in an event and try to match a track to them
    for (int cluster_n = 0; cluster_n < event->getNumberOfEcalClusters(); ++cluster_n) {
  
        // Get the cluster from the event 
        cluster = event->getEcalCluster(cluster_n);

        // Check that the cluster passes the energy requirement
        if (!passEnergyCut(cluster)) continue;

        // Check that the cluster passes the time requirement
        if (!passClusterTimeCut(cluster)) continue;

        // 
        if (!passClusterSizeCut(cluster)) continue;

        // Get the seed hit of the cluster
        EcalHit* seed_hit = cluster->getSeed();

        plotter->get1DHistogram("cluster energy")->Fill(cluster->getEnergy());
        plotter->get1DHistogram("cluster time")->Fill(cluster->getClusterTime());
        plotter->get2DHistogram("cluster energy v cluster time")->Fill(cluster->getEnergy(), 
                cluster->getClusterTime());
        plotter->get2DHistogram("cluster energy v cluster size")->Fill(cluster->getEnergy(), 
                cluster->getEcalHits()->GetEntriesFast());
        plotter->get2DHistogram("candidates")->Fill(
                seed_hit->getXCrystalIndex(), seed_hit->getYCrystalIndex(), 1);
        
        plotter->get1DHistogram("candidates - cluster energy")->Fill(cluster->getEnergy(), 1);
        plotter->get1DHistogram("candidates - cluster time")->Fill(cluster->getClusterTime(), 1);

        for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
            
            // Get a track from the event
            track = event->getTrack(track_n);

            plotter->get1DHistogram("track time")->Fill(track->getTrackTime());
            plotter->get1DHistogram("cluster time - track time")->Fill(
                    cluster->getClusterTime() - track->getTrackTime());

            if (isMatch(cluster, track)) {
                plotter->get2DHistogram("tracking efficiency")->Fill(seed_hit->getXCrystalIndex(), 
                    seed_hit->getYCrystalIndex(), 1);
                plotter->get1DHistogram("tracking efficiency - cluster energy")->Fill(cluster->getEnergy(), 1);
                plotter->get1DHistogram("tracking efficiency - cluster time")->Fill(cluster->getClusterTime(), 1);
                break; 
            } 
        }
    }
}

void SimpleTrackingEfficiencyAnalysis::finalize() { 
    plotter->get2DHistogram("tracking efficiency")->Divide(plotter->get2DHistogram("candidates"));
    plotter->get1DHistogram("tracking efficiency - cluster energy")->Divide(
            plotter->get1DHistogram("candidates - cluster energy"));
    plotter->get1DHistogram("tracking efficiency - cluster time")->Divide(
            plotter->get1DHistogram("candidates - cluster time"));

    plotter->saveToPdf("simple_tracking_efficiency.pdf");
}

void SimpleTrackingEfficiencyAnalysis::bookHistograms() { 

    plotter->build1DHistogram("cluster energy", 50, 0, 2);
    plotter->build1DHistogram("cluster time", 100, -20, 80);
    plotter->build1DHistogram("cluster time - track time", 100, -20, 80);
    plotter->build1DHistogram("track time", 100, -20, 80);
    plotter->build1DHistogram("candidates - cluster energy", 50, 0, 2);
    plotter->build1DHistogram("tracking efficiency - cluster energy", 50, 0, 2);
    plotter->build1DHistogram("candidates - cluster time", 100, -20, 80);
    plotter->build1DHistogram("tracking efficiency - cluster time", 100, -20, 80);

    plotter->build2DHistogram("cluster energy v cluster time", 50, 0, 2, 100, -20, 80);
    plotter->build2DHistogram("cluster energy v cluster size", 50, 0, 2, 10, 0, 10);
    plotter->build2DHistogram("candidates", 47, -23, 24, 12, -6, 6);
    plotter->build2DHistogram("tracking efficiency", 47, -23, 24, 12, -6, 6);
}

std::string SimpleTrackingEfficiencyAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;
}


bool SimpleTrackingEfficiencyAnalysis::passEnergyCut(EcalCluster* cluster) { 
    if (cluster->getEnergy() < cluster_energy_low_threshold
            || cluster->getEnergy() > cluster_energy_high_threshold) return false;

    return true;
}

bool SimpleTrackingEfficiencyAnalysis::passClusterTimeCut(EcalCluster* cluster) {   
    if (cluster->getClusterTime() < 40 || cluster->getClusterTime() > 50) return false;

    return true;   
}

bool SimpleTrackingEfficiencyAnalysis::passClusterSizeCut(EcalCluster* cluster) { 
    if (cluster->getEcalHits()->GetEntriesFast() < 3 
            || cluster->getEcalHits()->GetEntriesFast() > 6) return false;
    return true;

}

bool SimpleTrackingEfficiencyAnalysis::isMatch(EcalCluster* cluster, SvtTrack* track) { 

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
    if (abs(cluster_pos[0] - track_pos_at_cluster_shower_max[0]) > 20) return false;
    
    if (abs(cluster_pos[1] - track_pos_at_cluster_shower_max[1]) > 20) return false;

    return true;
}
