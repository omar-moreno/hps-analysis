
#include <TagProbeAnalysis.h>

TagProbeAnalysis::TagProbeAnalysis() 
    : plotter(new Plotter()), 
      total_events(0),
      total_trigger_events(0), 
      total_pair_events(0), 
      total_events_pass_fiducial_cut(0),
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
    //if (!event->isPair1Trigger()) return;
    //if (!event->isSingle1Trigger()) return;
    total_trigger_events++;

    // Search the event for a pair of clusters.  Only a time requirement is 
    // required at this point
    std::vector<EcalCluster*> pair = AnalysisUtils::getClusterPair(event);

    // If there aren't two clusters in the vector i.e. two matching clusters
    // weren't found, skip the event
    if (pair.size() != 2) return;
    total_pair_events++; 

    plotter->get1DHistogram("cluster time")->Fill(pair[0]->getClusterTime());
    plotter->get1DHistogram("cluster time")->Fill(pair[1]->getClusterTime());
    plotter->get2DHistogram("cluster pair energy")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    plotter->get2DHistogram("cluster pair time")->Fill(pair[0]->getClusterTime(), pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair dt")->Fill(pair[0]->getClusterTime() - pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair energy sum")->Fill(pair[0]->getEnergy() + pair[1]->getEnergy());

    plotter->get2DHistogram("cluster x position")->Fill(pair[0]->getPosition()[0], pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y position")->Fill(pair[0]->getPosition()[1], pair[1]->getPosition()[1]);

   
    // Check that the clusters aren't on the same side of the Ecal
    if (!passFiducialCut(pair[0], pair[1])) return;
    total_events_pass_fiducial_cut++;

    plotter->get1DHistogram("cluster time - cuts: fiducial")->Fill(pair[0]->getClusterTime());
    plotter->get1DHistogram("cluster time - cuts: fiducial")->Fill(pair[1]->getClusterTime());
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    plotter->get2DHistogram("cluster pair time - cuts: fiducial")->Fill(
            pair[0]->getClusterTime(), pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair dt - cuts: fiducial")->Fill(
            pair[0]->getClusterTime() - pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair energy sum - cuts: fiducial")->Fill(
            pair[0]->getEnergy() + pair[1]->getEnergy());

    plotter->get2DHistogram("cluster x position - cuts: fiducial")->Fill(pair[0]->getPosition()[0], pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y position - cuts: fiducial")->Fill(pair[0]->getPosition()[1], pair[1]->getPosition()[1]);


    if (!passClusterEnergySumCut(pair[0], pair[1])) return; 

    /*
    plotter->get1DHistogram("cluster time - cuts: fiducial, sum")->Fill(pair[0]->getClusterTime());
    plotter->get1DHistogram("cluster time - cuts: fiducial, sum")->Fill(pair[1]->getClusterTime());
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial, sum")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    plotter->get2DHistogram("cluster pair time - cuts: fiducial, sum")->Fill(
            pair[0]->getClusterTime(), pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair dt - cuts: fiducial, sum")->Fill(
            pair[0]->getClusterTime() - pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair energy sum - cuts: fiducial, sum")->Fill(
            pair[0]->getEnergy() + pair[1]->getEnergy());

    plotter->get2DHistogram("cluster x position - cuts: fiducial, sum")->Fill(pair[0]->getPosition()[0], pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y position - cuts: fiducial, sum")->Fill(pair[0]->getPosition()[1], pair[1]->getPosition()[1]);
    */

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
    
    // If the track is a positron, skip the event
    if (tag_track->getCharge() > 0) return;
   
    plotter->get1DHistogram("cluster time - candidates")->Fill(pair[0]->getClusterTime());
    plotter->get1DHistogram("cluster time - candidates")->Fill(pair[1]->getClusterTime());
    plotter->get2DHistogram("cluster pair energy - candidates")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    plotter->get2DHistogram("cluster pair time - candidates")->Fill(
            pair[0]->getClusterTime(), pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair dt - candidates")->Fill(
            pair[0]->getClusterTime() - pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair energy sum - candidates")->Fill(
            pair[0]->getEnergy() + pair[1]->getEnergy());

    total_tag_candidates++;

    plotter->get2DHistogram("tag clusters - candidates")->Fill( tag_seed_hit->getXCrystalIndex(), 
            tag_seed_hit->getYCrystalIndex(), 1); 

    plotter->get2DHistogram("cluster x position - candidates")->Fill(pair[0]->getPosition()[0], pair[1]->getPosition()[0]);

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
    if (probe_seed_hit->getYCrystalIndex() > 0) { 
        plotter->get1DHistogram("probe cluster energy - candidates - top")->Fill(probe_cluster->getEnergy());
    } else { 
        plotter->get1DHistogram("probe cluster energy - candidates - bottom")->Fill(probe_cluster->getEnergy());
    } 
        

    SvtTrack* probe_track = NULL; 
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
       
        // Don't attempt to match the the tag track
        if (tag_track == event->getTrack(track_n)) continue;

        // Get a track from the event
        probe_track = event->getTrack(track_n);
        
        if (this->isMatch(probe_cluster, probe_track)) { 
            plotter->get2DHistogram("probe clusters - matched")->Fill( probe_seed_hit->getXCrystalIndex(), 
                probe_seed_hit->getYCrystalIndex(), 1); 

            if (probe_seed_hit->getYCrystalIndex() > 0) { 
                plotter->get1DHistogram("probe cluster energy - matched - top")->Fill(probe_cluster->getEnergy());
            } else { 
                plotter->get1DHistogram("probe cluster energy - matched - bottom")->Fill(probe_cluster->getEnergy());
            } 

            double mass = AnalysisUtils::getInvariantMass(tag_track, probe_track);
            plotter->get1DHistogram("invariant mass - mollers")->Fill(mass); 
            break;
        }

        std::cout << "No match found" << std::endl;
        plotter->get2DHistogram("probe clusters - not matched")->Fill( probe_seed_hit->getXCrystalIndex(), 
                probe_seed_hit->getYCrystalIndex(), 1); 
        
        if (probe_seed_hit->getYCrystalIndex() > 0) { 
            plotter->get1DHistogram("probe cluster energy - not matched - top")->Fill(probe_cluster->getEnergy());
        } else { 
            plotter->get1DHistogram("probe cluster energy - not matched - bottom")->Fill(probe_cluster->getEnergy());
        }

        plotter->get2DHistogram("cluster x position - not matched")->Fill(pair[0]->getPosition()[0], pair[1]->getPosition()[0]);
        plotter->get2DHistogram("cluster y position - not matched")->Fill(pair[0]->getPosition()[1], pair[1]->getPosition()[1]);

        plotter->get2DHistogram("cluster pair energy - not matched")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    }
}

void TagProbeAnalysis::finalize() { 

    plotter->get2DHistogram("probe clusters - matched")->Divide(
            plotter->get2DHistogram("probe clusters - candidates"));

    plotter->setGraphType("asymm");
    ((TGraphAsymmErrors*) plotter->buildGraph("tracking efficiency - top"))->Divide(
        plotter->get1DHistogram("probe cluster energy - matched - top"),
        plotter->get1DHistogram("probe cluster energy - candidates - top"));

    ((TGraphAsymmErrors*) plotter->buildGraph("tracking efficiency - bottom"))->Divide(
        plotter->get1DHistogram("probe cluster energy - matched - bottom"),
        plotter->get1DHistogram("probe cluster energy - candidates - bottom"));

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

    plotter->build2DHistogram("cluster pair energy - cuts: fiducial", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial")->GetYaxis()->SetTitle("Cluster energy (GeV)");

    plotter->build2DHistogram("cluster pair energy - cuts: fiducial, sum", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial, sum")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial, sum")->GetYaxis()->SetTitle("Cluster energy (GeV)");

    plotter->build2DHistogram("cluster pair energy - candidates", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - candidates")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - candidates")->GetYaxis()->SetTitle("Cluster energy (GeV)");

    plotter->build2DHistogram("cluster pair energy - not matched", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - not matched")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - not matched")->GetYaxis()->SetTitle("Cluster energy (GeV)");

    plotter->build1DHistogram("Cluster pair energy sum", 50, 0, 1.5)->GetXaxis()->SetTitle(
            "Cluster pair energy sum (GeV)");
    plotter->build1DHistogram("Cluster pair energy sum - cuts: fiducial", 50, 0, 1.5)->GetXaxis()->SetTitle(
            "Cluster pair energy sum (GeV)");
    plotter->build1DHistogram("Cluster pair energy sum - cuts: fiducial, sum", 50, 0, 1.5)->GetXaxis()->SetTitle(
            "Cluster pair energy sum (GeV)");
    plotter->build1DHistogram("Cluster pair energy sum - candidates", 50, 0, 1.5)->GetXaxis()->SetTitle(
            "Cluster pair energy sum (GeV)");

    plotter->build1DHistogram("probe cluster energy - candidates - top", 50, 0, 1.5)->GetXaxis()->SetTitle("Probe Cluster Energy (GeV)");
    plotter->build1DHistogram("probe cluster energy - candidates - bottom", 50, 0, 1.5)->GetXaxis()->SetTitle("Probe Cluster Energy (GeV)");

    plotter->build1DHistogram("probe cluster energy - matched - top", 50, 0, 1.5)->GetXaxis()->SetTitle("Probe Cluster Energy (GeV)");
    plotter->build1DHistogram("probe cluster energy - matched - bottom", 50, 0, 1.5)->GetXaxis()->SetTitle("Probe Cluster Energy (GeV)");
    plotter->build1DHistogram("probe cluster energy - not matched - top", 50, 0, 1.5)->GetXaxis()->SetTitle("Probe Cluster Energy (GeV)");
    plotter->build1DHistogram("probe cluster energy - not matched - bottom", 50, 0, 1.5)->GetXaxis()->SetTitle("Probe Cluster Energy (GeV)");

    plotter->build2DHistogram("tag clusters", 47, -23, 24, 12, -6, 6);
    plotter->build2DHistogram("tag clusters - candidates", 47, -23, 24, 12, -6, 6);
    plotter->build2DHistogram("probe clusters - candidates", 47, -23, 24, 12, -6, 6);
    plotter->build2DHistogram("probe clusters - matched", 47, -23, 24, 12, -6, 6);
    plotter->build2DHistogram("probe clusters - not matched", 47, -23, 24, 12, -6, 6);



    // Cluster time //
    
    plotter->build1DHistogram("cluster time", 160, 20, 100)->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->build1DHistogram("cluster time - cuts: fiducial", 160, 20, 100)->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->build1DHistogram("cluster time - cuts: fiducial, sum", 160, 20, 100)->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->build1DHistogram("cluster time - candidates", 160, 20, 100)->GetXaxis()->SetTitle("Cluster time (ns)");

    plotter->build2DHistogram("cluster pair time", 160, 20, 100, 160, 20, 100);
    plotter->get2DHistogram("cluster pair time")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build2DHistogram("cluster pair time - cuts: fiducial", 160, 20, 100, 160, 20, 100);
    plotter->get2DHistogram("cluster pair time - cuts: fiducial")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - cuts: fiducial")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build2DHistogram("cluster pair time - cuts: fiducial, sum", 160, 20, 100, 160, 20, 100);
    plotter->get2DHistogram("cluster pair time - cuts: fiducial, sum")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - cuts: fiducial, sum")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build2DHistogram("cluster pair time - candidates", 160, 20, 100, 160, 20, 100);
    plotter->get2DHistogram("cluster pair time - candidates")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - candidates")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build1DHistogram("Cluster pair dt", 100, -5, 5)->GetXaxis()->SetTitle(
            "Cluster pair dt");

    plotter->build1DHistogram("Cluster pair dt - cuts: fiducial", 100, -5, 5)->GetXaxis()->SetTitle(
            "Cluster pair dt");
    
    plotter->build1DHistogram("Cluster pair dt - cuts: fiducial, sum", 100, -5, 5)->GetXaxis()->SetTitle(
            "Cluster pair dt");

    plotter->build1DHistogram("Cluster pair dt - candidates", 100, -5, 5)->GetXaxis()->SetTitle(
            "Cluster pair dt");

    // Cluster position
    plotter->build2DHistogram("cluster x position", 100, -200, 100, 100, -200, 100);
    plotter->get2DHistogram("cluster x position")->GetXaxis()->SetTitle("First cluster x position (mm)");
    plotter->get2DHistogram("cluster x position")->GetYaxis()->SetTitle("Second cluster x position (mm)");
    
    plotter->build2DHistogram("cluster y position", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster y position")->GetXaxis()->SetTitle("First cluster y position (mm)");
    plotter->get2DHistogram("cluster y position")->GetYaxis()->SetTitle("Second cluster y position (mm)");

    plotter->build2DHistogram("cluster x position - cuts: fiducial", 100, -200, 100, 100, -200, 100);
    plotter->get2DHistogram("cluster x position - cuts: fiducial")->GetXaxis()->SetTitle("First cluster x position (mm)");
    plotter->get2DHistogram("cluster x position - cuts: fiducial")->GetYaxis()->SetTitle("Second cluster x position (mm)");

    plotter->build2DHistogram("cluster x position - candidates", 100, -200, 100, 100, -200, 100);
    plotter->get2DHistogram("cluster x position - candidates")->GetXaxis()->SetTitle("First cluster x position (mm)");
    plotter->get2DHistogram("cluster x position - candidates")->GetYaxis()->SetTitle("Second cluster x position (mm)");

    plotter->build2DHistogram("cluster y position - cuts: fiducial", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster y position - cuts: fiducial")->GetXaxis()->SetTitle("First cluster y position (mm)");
    plotter->get2DHistogram("cluster y position - cuts: fiducial")->GetYaxis()->SetTitle("Second cluster y position (mm)");

    plotter->build2DHistogram("cluster x position - cuts: fiducial, sum", 100, -200, 100, 100, -200, 100);
    plotter->get2DHistogram("cluster x position - cuts: fiducial, sum")->GetXaxis()->SetTitle("First cluster x position (mm)");
    plotter->get2DHistogram("cluster x position - cuts: fiducial, sum")->GetYaxis()->SetTitle("Second cluster x position (mm)");
    
    plotter->build2DHistogram("cluster y position - cuts: fiducial, sum", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster y position - cuts: fiducial, sum")->GetXaxis()->SetTitle("First cluster y position (mm)");
    plotter->get2DHistogram("cluster y position - cuts: fiducial, sum")->GetYaxis()->SetTitle("Second cluster y position (mm)");

    plotter->build2DHistogram("cluster x position - not matched", 100, -200, 100, 100, -200, 100);
    plotter->get2DHistogram("cluster x position - not matched")->GetXaxis()->SetTitle("First cluster x position (mm)");
    plotter->get2DHistogram("cluster x position - not matched")->GetYaxis()->SetTitle("Second cluster x position (mm)");
    
    plotter->build2DHistogram("cluster y position - not matched", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster y position - not matched")->GetXaxis()->SetTitle("First cluster y position (mm)");
    plotter->get2DHistogram("cluster y position - not matched")->GetYaxis()->SetTitle("Second cluster y position (mm)");

    // Invariant mass
    plotter->build1DHistogram("invariant mass - mollers", 50, 0, 0.1)->GetXaxis()->SetTitle("Invariant mass (GeV)");

}

std::string TagProbeAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;
}

bool TagProbeAnalysis::passClusterEnergySumCut(EcalCluster* first_cluster, EcalCluster* second_cluster) {
    
    double cluster_energy_sum = first_cluster->getEnergy() + second_cluster->getEnergy(); 
    if (cluster_energy_sum > 1.1 || cluster_energy_sum < .85) return false; 

    return true;
}

bool TagProbeAnalysis::passFiducialCut(EcalCluster* first_cluster, EcalCluster* second_cluster) { 
   
    // Make sure that the clusters aren't on the same side in y
    if (first_cluster->getPosition()[1]*second_cluster->getPosition()[1] > 0) return false;

    // Require that they are on the electron side in x
    if (first_cluster->getPosition()[0] > 0 || second_cluster->getPosition()[0] > 0) return false;

    if (first_cluster->getPosition()[0]+second_cluster->getPosition()[0] > -135) return false;
    
    if (first_cluster->getPosition()[0]+second_cluster->getPosition()[0] < -190) return false;

    return true;
}


bool TagProbeAnalysis::isMatch(EcalCluster* cluster, SvtTrack* track) { 

    // Check that the track and cluster are in the same detector volume.
    // If not, thre is no way they can match.
    if (track->isTopTrack() && cluster->getPosition()[1] < 0
            || track->isBottomTrack() && cluster->getPosition()[1] > 0) return false;

    // Get the cluster position
    std::vector<double> cluster_pos = cluster->getPosition();
    /*std::cout << "[ TagProbeAnalysis ]: ECal cluster position: " 
        << " x: " << cluster_pos[0] 
        << " y: " << cluster_pos[1] 
        << " z: " << cluster_pos[2]
        << std::endl;*/

    // Extrapolate the track to the Ecal cluster position 
    std::vector<double> track_pos_at_ecal 
        = TrackExtrapolator::extrapolateTrack(track, cluster_pos[2]);
    /*std::cout << "[ TagProbeAnalysis ]: Track position at shower max: " 
         << " x: " << track_pos_at_ecal[0] 
         << " y: " << track_pos_at_ecal[1] 
         << " z: " << track_pos_at_ecal[2]
         << std::endl;*/ 

    double p = AnalysisUtils::getMagnitude(track->getMomentum());

    double delta_x = cluster_pos[0] - track_pos_at_ecal[0];
    double delta_y = cluster_pos[1] - track_pos_at_ecal[1];

    // Check that dx and dy between the extrapolated track and cluster
    // positions is reasonable
    if ((track->isTopTrack() && (delta_x > 14 || delta_x < -18)) ||
        (track->isBottomTrack() && (delta_x > 9 || delta_x < -21))) return false;

    if ((track->isTopTrack() && (delta_y > 14 || delta_y < -14)) ||
        (track->isBottomTrack() && (delta_y > 14 || delta_y < -14))) return false;

    return true;
}
