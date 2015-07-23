
#include <TagProbeAnalysis.h>

TagProbeAnalysis::TagProbeAnalysis() 
    : plotter(new Plotter()), 
      total_events(0),
      total_pair_trigger_events(0), 
      total_pair_events(0), 
      total_tag_candidates(0), 
      class_name("TagProbeAnalysis") {       
}

TagProbeAnalysis::~TagProbeAnalysis() { 
}

void TagProbeAnalysis::initialize() { 
    this->bookHistograms();
    
    srand(time(0)); 
}

void TagProbeAnalysis::processEvent(HpsEvent* event) {

    total_events++;

    // Only look at pair1 triggers
    if (!event->isPair1Trigger()) return;
    total_pair_trigger_events++;

    // Search the event for a pair of clusters.  Only a time requirement is 
    // required at this point
    std::vector<EcalCluster*> pair = AnalysisUtils::getClusterPair(event);

    // If there aren't two clusters in the vector i.e. two matching clusters
    // weren't found, skip the event
    if (pair.size() != 2) return;
    total_pair_events++; 

    plotter->get2DHistogram("cluster pair energy")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    plotter->get2DHistogram("cluster pair time")->Fill(pair[0]->getClusterTime(), pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair dt")->Fill(pair[0]->getClusterTime() - pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair energy sum")->Fill(pair[0]->getEnergy() + pair[1]->getEnergy());

    plotter->get2DHistogram("cluster x position")->Fill(pair[0]->getPosition()[0], pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y position")->Fill(pair[0]->getPosition()[1], pair[1]->getPosition()[1]);

   
    // Check that the clusters aren't on the same side of the Ecal
    if (!passFiducialCut(pair[0], pair[1])) return;

    plotter->get2DHistogram("cluster pair energy - cuts: fiducial")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    plotter->get2DHistogram("cluster pair time - cuts: fiducial")->Fill(
            pair[0]->getClusterTime(), pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair dt - cuts: fiducial")->Fill(
            pair[0]->getClusterTime() - pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair energy sum - cuts: fiducial")->Fill(
            pair[0]->getEnergy() + pair[1]->getEnergy());

    plotter->get2DHistogram("cluster x position - cuts: fiducial")->Fill(pair[0]->getPosition()[0], pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y position - cuts: fiducial")->Fill(pair[0]->getPosition()[1], pair[1]->getPosition()[1]);


    // Check if the clusters pass the cluster energy sum cut
    if (!passClusterEnergySumCut(pair[0], pair[1])) return; 

    plotter->get2DHistogram("cluster pair energy - cuts: fiducial, sum")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    plotter->get2DHistogram("cluster pair time - cuts: fiducial, sum")->Fill(
            pair[0]->getClusterTime(), pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair dt - cuts: fiducial, sum")->Fill(
            pair[0]->getClusterTime() - pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair energy sum - cuts: fiducial, sum")->Fill(
            pair[0]->getEnergy() + pair[1]->getEnergy());

    plotter->get2DHistogram("cluster x position - cuts: fiducial, sum")->Fill(pair[0]->getPosition()[0], pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y position - cuts: fiducial, sum")->Fill(pair[0]->getPosition()[1], pair[1]->getPosition()[1]);

    // Randomly choose one of the two ECal clusters
    double cluster_index = rand()%2;
    EcalCluster* tag_cluster = pair[cluster_index]; 

    EcalHit* tag_seed_hit = tag_cluster->getSeed();

    plotter->get2DHistogram("tag clusters")->Fill( tag_seed_hit->getXCrystalIndex(), 
            tag_seed_hit->getYCrystalIndex(), 1); 

    // Try and match a track to the cluster.
    SvtTrack* tag_track = NULL;
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
        
        if (this->isMatch(tag_cluster, event->getTrack(track_n))) { 
           tag_track = event->getTrack(track_n); 
           break;
        }
    }

    // If no match was found, skip the event
    if (tag_track == NULL) return;
    total_tag_candidates++;

    plotter->get2DHistogram("tag clusters - candidates")->Fill( tag_seed_hit->getXCrystalIndex(), 
            tag_seed_hit->getYCrystalIndex(), 1); 

    // Get the probe cluster
    EcalCluster* probe_cluster = NULL;
    if (cluster_index == 0) { 
        probe_cluster = pair[1];
    } else {
        probe_cluster = pair[0];
    }

    EcalHit* probe_seed_hit = probe_cluster->getSeed();
    plotter->get2DHistogram("probe clusters - candidates")->Fill( probe_seed_hit->getXCrystalIndex(), 
            probe_seed_hit->getYCrystalIndex(), 1); 
   
    SvtTrack* probe_track = NULL; 
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
       
        // Don't attempt to match the the tag track
        if (tag_track == event->getTrack(track_n)) continue;

        // Get a track from the event
        probe_track = event->getTrack(track_n);
        
        if (this->isMatch(probe_cluster, probe_track)) { 
            plotter->get2DHistogram("probe clusters - matched")->Fill( probe_seed_hit->getXCrystalIndex(), 
                probe_seed_hit->getYCrystalIndex(), 1); 
            break;
        }
    } 
}

void TagProbeAnalysis::finalize() { 

    plotter->get2DHistogram("probe clusters - matched")->Divide(
            plotter->get2DHistogram("probe clusters - candidates"));

    plotter->saveToPdf("tag_probe_efficiency.pdf");
    plotter->saveToRootFile("tag_probe_efficiency.root");

    std::cout << "//---------------------------------------------------//" << std::endl;
    std::cout << "// Total events: " << total_events << std::endl;
    std::cout << "//---------------------------------------------------//" << std::endl;

    return;
}

void TagProbeAnalysis::bookHistograms() { 
 
    plotter->setType("float");

    // Cluster energy // 
    plotter->build2DHistogram("cluster pair energy", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy")->GetYaxis()->SetTitle("Cluster energy (GeV)");

    plotter->build2DHistogram("cluster pair energy - cuts: fiducial, sum", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial, sum")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial, sum")->GetYaxis()->SetTitle("Cluster energy (GeV)");

    plotter->build2DHistogram("cluster pair energy - cuts: fiducial", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial")->GetYaxis()->SetTitle("Cluster energy (GeV)");


    plotter->build1DHistogram("Cluster pair energy sum", 50, 0, 1.5)->GetXaxis()->SetTitle(
            "Cluster pair energy sum (GeV)");
    plotter->build1DHistogram("Cluster pair energy sum - cuts: fiducial, sum", 50, 0, 1.5)->GetXaxis()->SetTitle(
            "Cluster pair energy sum (GeV)");
    plotter->build1DHistogram("Cluster pair energy sum - cuts: fiducial", 50, 0, 1.5)->GetXaxis()->SetTitle(
            "Cluster pair energy sum (GeV)");


    plotter->build2DHistogram("tag clusters", 47, -23, 24, 12, -6, 6);
    plotter->build2DHistogram("tag clusters - candidates", 47, -23, 24, 12, -6, 6);
    plotter->build2DHistogram("probe clusters - candidates", 47, -23, 24, 12, -6, 6);
    plotter->build2DHistogram("probe clusters - matched", 47, -23, 24, 12, -6, 6);

    // Cluster time //
    plotter->build2DHistogram("cluster pair time", 160, 20, 100, 160, 20, 100);
    plotter->get2DHistogram("cluster pair time")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build2DHistogram("cluster pair time - cuts: fiducial, sum", 160, 20, 100, 160, 20, 100);
    plotter->get2DHistogram("cluster pair time - cuts: fiducial, sum")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - cuts: fiducial, sum")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build2DHistogram("cluster pair time - cuts: fiducial", 160, 20, 100, 160, 20, 100);
    plotter->get2DHistogram("cluster pair time - cuts: fiducial")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - cuts: fiducial")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build1DHistogram("Cluster pair dt", 100, -5, 5)->GetXaxis()->SetTitle(
            "Cluster pair dt");

    plotter->build1DHistogram("Cluster pair dt - cuts: fiducial, sum", 100, -5, 5)->GetXaxis()->SetTitle(
            "Cluster pair dt");

    plotter->build1DHistogram("Cluster pair dt - cuts: fiducial", 100, -5, 5)->GetXaxis()->SetTitle(
            "Cluster pair dt");

    // Cluster position
    plotter->build2DHistogram("cluster x position", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster x position")->GetXaxis()->SetTitle("First cluster x position (mm)");
    plotter->get2DHistogram("cluster x position")->GetYaxis()->SetTitle("Second cluster x position (mm)");
    
    plotter->build2DHistogram("cluster y position", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster y position")->GetXaxis()->SetTitle("First cluster y position (mm)");
    plotter->get2DHistogram("cluster y position")->GetYaxis()->SetTitle("Second cluster y position (mm)");

    plotter->build2DHistogram("cluster x position - cuts: fiducial", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster x position - cuts: fiducial")->GetXaxis()->SetTitle("First cluster x position (mm)");
    plotter->get2DHistogram("cluster x position - cuts: fiducial")->GetYaxis()->SetTitle("Second cluster x position (mm)");
    
    plotter->build2DHistogram("cluster y position - cuts: fiducial", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster y position - cuts: fiducial")->GetXaxis()->SetTitle("First cluster y position (mm)");
    plotter->get2DHistogram("cluster y position - cuts: fiducial")->GetYaxis()->SetTitle("Second cluster y position (mm)");

    plotter->build2DHistogram("cluster x position - cuts: fiducial, sum", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster x position - cuts: fiducial, sum")->GetXaxis()->SetTitle("First cluster x position (mm)");
    plotter->get2DHistogram("cluster x position - cuts: fiducial, sum")->GetYaxis()->SetTitle("Second cluster x position (mm)");
    
    plotter->build2DHistogram("cluster y position - cuts: fiducial, sum", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster y position - cuts: fiducial, sum")->GetXaxis()->SetTitle("First cluster y position (mm)");
    plotter->get2DHistogram("cluster y position - cuts: fiducial, sum")->GetYaxis()->SetTitle("Second cluster y position (mm)");

    /*
    ecal_plotter->setType("float");
    ecal_plotter->build1DHistogram("Cluster pair energy sum - Pass Cuts", 50, 0, 2);
    ecal_plotter->build1DHistogram("Cluster pair dt - Pass cuts", 100, -50, 50);
    ecal_plotter->build1DHistogram("Cluster index", 10, 0, 10);
    ecal_plotter->build1DHistogram("Tag cluster time - track time", 100, -20, 80);
    ecal_plotter->build1DHistogram("Tag cluster time - track time - Pass cuts", 100, -20, 80);
    ecal_plotter->build1DHistogram("Probe cluster time - track time", 100, -20, 80);
    ecal_plotter->build1DHistogram("Probe cluster time - track time - Pass cuts", 100, -20, 80);
    ecal_plotter->build1DHistogram("E/p", 100, 0, 2);
    ecal_plotter->build1DHistogram("E/p - Tag", 100, 0, 2);
    ecal_plotter->build1DHistogram("E/p - Probe", 100, 0, 2);
    ecal_plotter->build1DHistogram("Tag cluster energy", 50, 0, 1);
    ecal_plotter->build2DHistogram("Tag cluster energy vs Cluster x index", 47, -23, 24, 50, 0, 1);
    ecal_plotter->build2DHistogram("Tag cluster energy vs Cluster y index", 50, 0, 1, 12, -6, 6);
    ecal_plotter->build1DHistogram("Probe cluster energy", 50, 0, 1);
    ecal_plotter->build2DHistogram("Probe cluster energy vs Cluster x index", 47, -23, 24, 50, 0, 1);
    ecal_plotter->build2DHistogram("Probe cluster energy vs Cluster y index", 50, 0, 1, 12, -6, 6);
    ecal_plotter->build2DHistogram("Probe clusters", 47, -23, 24, 12, -6, 6);
    ecal_plotter->build2DHistogram("Matched clusters", 47, -23, 24, 12, -6, 6);
    ecal_plotter->build2DHistogram("Tracking efficiency", 47, -23, 24, 12, -6, 6);

    svt_plotter->setType("float");
    svt_plotter->build1DHistogram("Track-Cluster dx", 100, -200, 200);
    svt_plotter->build1DHistogram("Track-Cluster dy", 100, -200, 200);
    svt_plotter->build1DHistogram("Track-Cluster dx - Pass cuts", 100, -200, 200);
    svt_plotter->build1DHistogram("Track-Cluster dy - Pass cuts", 100, -200, 200);
    */
    return;
}

std::string TagProbeAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;
}

bool TagProbeAnalysis::passClusterEnergySumCut(EcalCluster* first_cluster, EcalCluster* second_cluster) {
    
    double cluster_energy_sum = first_cluster->getEnergy() + second_cluster->getEnergy(); 
    if (cluster_energy_sum > 1.15 || cluster_energy_sum < .7) return false; 

    return true;
}

bool TagProbeAnalysis::passFiducialCut(EcalCluster* first_cluster, EcalCluster* second_cluster) { 
   
    // Make sure that the clusters aren't on the same side in y
    if (first_cluster->getPosition()[1]*second_cluster->getPosition()[1] > 0) return false;

    // Require that they are on the electron side in x
    if (first_cluster->getPosition()[0] > 0 && second_cluster->getPosition()[0] > 0) return false;

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


    double p = AnalysisUtils::getMagnitude(track->getMomentum());

    // If the track and cluster are in opposite volumes, then they can't 
    // be a match
    if (cluster_pos[1]*track_pos_at_cluster_shower_max[1] < 0) return false;

    if (track->isTopTrack()) { 
        /*plotter->get2DHistogram("cluster x v extrapolated track x - top")->Fill(cluster_pos[0], 
                track_pos_at_cluster_shower_max[0]);
        plotter->get2DHistogram("cluster y v extrapolated track y - top")->Fill(cluster_pos[1], 
                track_pos_at_cluster_shower_max[1]);

        plotter->get1DHistogram("cluster x - extrapolated track x - top")->Fill(cluster_pos[0] 
                - track_pos_at_cluster_shower_max[0]);
        plotter->get1DHistogram("cluster y - extrapolated track y - top")->Fill(cluster_pos[1] 
                - track_pos_at_cluster_shower_max[1]);

        plotter->get2DHistogram("p v extrapolated track x - top")->Fill(p, track_pos_at_cluster_shower_max[0]);
        plotter->get2DHistogram("p v cluster x - top")->Fill(p, cluster_pos[0]);
        plotter->get2DHistogram("cluster pair energy v cluster x - top")->Fill(cluster->getEnergy(), cluster_pos[0]);
        plotter->get2DHistogram("cluster x - track x v e/p - top")->Fill(cluster_pos[0] - track_pos_at_cluster_shower_max[0],
               cluster->getEnergy()/p); 
        plotter->get2DHistogram("cluster x v e/p - top")->Fill(cluster_pos[0],
               cluster->getEnergy()/p); 
        plotter->get2DHistogram("cluster y v e/p - top")->Fill(cluster_pos[1],
               cluster->getEnergy()/p); 
        */
        if (track->getCharge() < 0) { 
          //  plotter->get2DHistogram("cluster x v extrapolated track x - top - electrons")->Fill(cluster_pos[0], 
          //          track_pos_at_cluster_shower_max[0]);
        } else {
          //  plotter->get2DHistogram("cluster x v extrapolated track x - top - positrons")->Fill(cluster_pos[0], 
          //          track_pos_at_cluster_shower_max[0]);
        }
    
    } else {
       /* plotter->get2DHistogram("cluster x v extrapolated track x - bottom")->Fill(cluster_pos[0], 
                track_pos_at_cluster_shower_max[0]);
        plotter->get2DHistogram("cluster y v extrapolated track y - bottom")->Fill(cluster_pos[1], 
                track_pos_at_cluster_shower_max[1]);

        plotter->get1DHistogram("cluster x - extrapolated track x - bottom")->Fill(cluster_pos[0] 
                - track_pos_at_cluster_shower_max[0]);
        plotter->get1DHistogram("cluster y - extrapolated track y - bottom")->Fill(cluster_pos[1] 
                - track_pos_at_cluster_shower_max[1]);
        
        plotter->get2DHistogram("p v extrapolated track x - bottom")->Fill(p, track_pos_at_cluster_shower_max[0]);
        plotter->get2DHistogram("p v cluster x - bottom")->Fill(p, cluster_pos[0]);
        plotter->get2DHistogram("cluster pair energy v cluster x - bottom")->Fill(cluster->getEnergy(), cluster_pos[0]);
        plotter->get2DHistogram("cluster x - track x v e/p - bottom")->Fill(cluster_pos[0] - track_pos_at_cluster_shower_max[0],
               cluster->getEnergy()/p); 
        plotter->get2DHistogram("cluster x v e/p - bottom")->Fill(cluster_pos[0],
               cluster->getEnergy()/p); 
        plotter->get2DHistogram("cluster y v e/p - bottom")->Fill(cluster_pos[1],
               cluster->getEnergy()/p); 
        
        if (track->getCharge() < 0) { 
            plotter->get2DHistogram("cluster x v extrapolated track x - bottom - electrons")->Fill(cluster_pos[0], 
                    track_pos_at_cluster_shower_max[0]);
        } else {
            plotter->get2DHistogram("cluster x v extrapolated track x - bottom - positrons")->Fill(cluster_pos[0], 
                    track_pos_at_cluster_shower_max[0]);
        }*/
    }
    
    // Check that dx and dy between the extrapolated track and cluster
    // positions is reasonable
    if (std::abs(cluster_pos[0] - track_pos_at_cluster_shower_max[0]) > 20) return false;

    if (cluster->getEnergy()/p < .5) return false;
    
    /*
    if (track->isTopTrack()) { 
        plotter->get2DHistogram("cluster x v extrapolated track x - top - matched")->Fill(cluster_pos[0], 
                track_pos_at_cluster_shower_max[0]);
        plotter->get2DHistogram("cluster y v extrapolated track y - top - matched")->Fill(cluster_pos[1], 
                track_pos_at_cluster_shower_max[1]);

        plotter->get1DHistogram("cluster x - extrapolated track x - top - matched")->Fill(cluster_pos[0] 
                - track_pos_at_cluster_shower_max[0]);
        plotter->get1DHistogram("cluster y - extrapolated track y - top - matched")->Fill(cluster_pos[1] 
                - track_pos_at_cluster_shower_max[1]);
    } else {
        plotter->get2DHistogram("cluster x v extrapolated track x - bottom - matched")->Fill(cluster_pos[0], 
                track_pos_at_cluster_shower_max[0]);
        plotter->get2DHistogram("cluster y v extrapolated track y - bottom - matched")->Fill(cluster_pos[1], 
                track_pos_at_cluster_shower_max[1]);

        plotter->get1DHistogram("cluster x - extrapolated track x - bottom - matched")->Fill(cluster_pos[0] 
                - track_pos_at_cluster_shower_max[0]);
        plotter->get1DHistogram("cluster y - extrapolated track y - bottom - matched")->Fill(cluster_pos[1] 
                - track_pos_at_cluster_shower_max[1]);
    }*/

    /*if (cluster_pos[0] - track_pos_at_cluster_shower_max[0] > 30 ||
            cluster_pos[0] - track_pos_at_cluster_shower_max[0] < -30) return false;
    
    if (cluster_pos[1] - track_pos_at_cluster_shower_max[1] > 30
            || cluster_pos[1] - track_pos_at_cluster_shower_max[1] < -30) return false;*/

    //std::cout << "Track and cluster are a match" << std::endl;

    return true;
}
