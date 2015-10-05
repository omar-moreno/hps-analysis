
#include <TagProbeAnalysis.h>

TagProbeAnalysis::TagProbeAnalysis() 
    : plotter(new Plotter()),
      ecal_utils(new EcalUtils()), 
      matcher(new TrackClusterMatcher()),  
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


    // Only look at single 1 triggers
    /*if (!event->isSingle1Trigger()) return;

    // Only look at events with the SVT bias ON
    if (!event->isSvtBiasOn()) return; 
    
    // Only look at events where the SVT is closed
    if (!event->isSvtClosed()) return;
    */

    // Get a "good" pair from the event.  If a good pair isn't found, skip
    // the event.
    std::vector<EcalCluster*> pair = ecal_utils->getClusterPair(event);
    if (pair.size() != 2 || pair[0] == nullptr || pair[1] == nullptr) return;

    total_pair_events++; 

    // Check that the clusters aren't on the same side of the Ecal
    if (!passFiducialCut(pair[0], pair[1])) return;
    total_events_pass_fiducial_cut++;

    plotter->get2DHistogram("cluster pair energy - cuts: fiducial")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    plotter->get2DHistogram("cluster pair time - cuts: fiducial")->Fill(
            pair[0]->getClusterTime(), pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair dt - cuts: fiducial")->Fill(
            pair[0]->getClusterTime() - pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair energy sum - cuts: fiducial")->Fill(
            pair[0]->getEnergy() + pair[1]->getEnergy());
    plotter->get2DHistogram("cluster x vs cluster x - cuts: fiducial")->Fill(
            pair[0]->getPosition()[0], pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y vs cluster y - cuts: fiducial")->Fill(
            pair[0]->getPosition()[1], pair[1]->getPosition()[1]);

    if (!passClusterEnergySumCut(pair[0], pair[1])) return; 
 
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial, sum")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    plotter->get2DHistogram("cluster pair time - cuts: fiducial, sum")->Fill(
            pair[0]->getClusterTime(), pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair dt - cuts: fiducial, sum")->Fill(
            pair[0]->getClusterTime() - pair[1]->getClusterTime());
    plotter->get1DHistogram("Cluster pair energy sum - cuts: fiducial, sum")->Fill(
            pair[0]->getEnergy() + pair[1]->getEnergy());
    plotter->get2DHistogram("cluster x vs cluster x - cuts: fiducial, sum")->Fill(
            pair[0]->getPosition()[0], pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y vs cluster y - cuts: fiducial, sum")->Fill(
            pair[0]->getPosition()[1], pair[1]->getPosition()[1]);

    // Randomly choose one of the two ECal clusters
    double cluster_index = rand()%2;
    EcalCluster* tag_cluster = pair[cluster_index]; 

    EcalHit* tag_seed_hit = tag_cluster->getSeed();

    plotter->get2DHistogram("tag clusters")->Fill( tag_seed_hit->getXCrystalIndex(), 
            tag_seed_hit->getYCrystalIndex(), 1); 

    // Find all track-cluster matches in the event
    matcher->findAllMatches(event);

    // Check if matches was found for the two clusters
    SvtTrack* tag_track = matcher->getMatchingTrack(tag_cluster); 
    if (tag_track == nullptr) return; 

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

    plotter->get2DHistogram("cluster x vs cluster x - candidates")->Fill(pair[0]->getPosition()[0], pair[1]->getPosition()[0]);

    double p = AnalysisUtils::getMagnitude(tag_track->getMomentum());
    plotter->get2DHistogram("p v theta - candidates")->Fill(p, 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(tag_track))));


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

    // Check if matches was found for the two clusters
    SvtTrack* probe_track = matcher->getMatchingTrack(probe_cluster); 
    if (probe_track == nullptr) { 
        
        plotter->get2DHistogram("probe clusters - not matched")->Fill( probe_seed_hit->getXCrystalIndex(), 
                probe_seed_hit->getYCrystalIndex(), 1); 
        
        if (probe_seed_hit->getYCrystalIndex() > 0) { 
            plotter->get1DHistogram("probe cluster energy - not matched - top")->Fill(probe_cluster->getEnergy());
        } else { 
            plotter->get1DHistogram("probe cluster energy - not matched - bottom")->Fill(probe_cluster->getEnergy());
        }
        plotter->get2DHistogram("cluster x vs cluster x - not matched")->Fill(pair[0]->getPosition()[0],
                pair[1]->getPosition()[0]);
        plotter->get2DHistogram("cluster y position - not matched")->Fill(pair[0]->getPosition()[1],
                pair[1]->getPosition()[1]);
        plotter->get2DHistogram("cluster pair energy - not matched")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());

        p = AnalysisUtils::getMagnitude(tag_track->getMomentum());
        plotter->get2DHistogram("p v theta - not matched")->Fill(p, 
                std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(tag_track))));


        return;
    } 

    plotter->get2DHistogram("probe clusters - matched")->Fill( probe_seed_hit->getXCrystalIndex(), 
            probe_seed_hit->getYCrystalIndex(), 1); 
            
    if (probe_seed_hit->getYCrystalIndex() > 0) { 
        plotter->get1DHistogram("probe cluster energy - matched - top")->Fill(probe_cluster->getEnergy());
    } else { 
        plotter->get1DHistogram("probe cluster energy - matched - bottom")->Fill(probe_cluster->getEnergy());
    } 

    p = AnalysisUtils::getMagnitude(tag_track->getMomentum());
    plotter->get2DHistogram("p v theta - matched")->Fill(p, 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(tag_track))));
    

    double mass = AnalysisUtils::getInvariantMass(tag_track, probe_track);
    plotter->get1DHistogram("invariant mass - mollers")->Fill(mass); 
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

    ecal_utils->saveHistograms();
    
    std::cout << "//---------------------------------------------------//" << std::endl;
    std::cout << "// Total events: " << total_events << std::endl;
    std::cout << "//---------------------------------------------------//" << std::endl;

    return;
}

void TagProbeAnalysis::bookHistograms() { 
 
    plotter->setType("float");

    // Cluster energy // 
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
    
    plotter->build1DHistogram("cluster time - cuts: fiducial, sum", 160, 20, 100)->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->build1DHistogram("cluster time - candidates", 160, 20, 100)->GetXaxis()->SetTitle("Cluster time (ns)");

    plotter->build2DHistogram("cluster pair time - cuts: fiducial", 160, 20, 100, 160, 20, 100);
    plotter->get2DHistogram("cluster pair time - cuts: fiducial")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - cuts: fiducial")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build2DHistogram("cluster pair time - cuts: fiducial, sum", 160, 20, 100, 160, 20, 100);
    plotter->get2DHistogram("cluster pair time - cuts: fiducial, sum")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - cuts: fiducial, sum")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build2DHistogram("cluster pair time - candidates", 160, 20, 100, 160, 20, 100);
    plotter->get2DHistogram("cluster pair time - candidates")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - candidates")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build1DHistogram("Cluster pair dt - cuts: fiducial", 100, -5, 5)->GetXaxis()->SetTitle(
            "Cluster pair dt");
    
    plotter->build1DHistogram("Cluster pair dt - cuts: fiducial, sum", 100, -5, 5)->GetXaxis()->SetTitle(
            "Cluster pair dt");

    plotter->build1DHistogram("Cluster pair dt - candidates", 100, -5, 5)->GetXaxis()->SetTitle(
            "Cluster pair dt");

    // Cluster position
    plotter->build2DHistogram("cluster x vs cluster x - cuts: fiducial", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster x vs cluster x - cuts: fiducial")->GetXaxis()->SetTitle("Cluster position - x (mm)");
    plotter->get2DHistogram("cluster x vs cluster x - cuts: fiducial")->GetYaxis()->SetTitle("Cluster position - x (mm)");
    
    plotter->build2DHistogram("cluster x vs cluster x - candidates", 100, -200, 100, 100, -200, 100);
    plotter->get2DHistogram("cluster x vs cluster x - candidates")->GetXaxis()->SetTitle("First cluster x vs cluster x (mm)");
    plotter->get2DHistogram("cluster x vs cluster x - candidates")->GetYaxis()->SetTitle("Second cluster x vs cluster x (mm)");

    plotter->build2DHistogram("cluster y vs cluster y - cuts: fiducial", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster y vs cluster y - cuts: fiducial")->GetXaxis()->SetTitle("First cluster y position (mm)");
    plotter->get2DHistogram("cluster y vs cluster y - cuts: fiducial")->GetYaxis()->SetTitle("Second cluster y position (mm)");

    plotter->build2DHistogram("cluster x vs cluster x - cuts: fiducial, sum", 100, -200, 100, 100, -200, 100);
    plotter->get2DHistogram("cluster x vs cluster x - cuts: fiducial, sum")->GetXaxis()->SetTitle("First cluster x position (mm)");
    plotter->get2DHistogram("cluster x vs cluster x - cuts: fiducial, sum")->GetYaxis()->SetTitle("Second cluster x position (mm)");
    
    plotter->build2DHistogram("cluster y vs cluster y - cuts: fiducial, sum", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster y vs cluster y - cuts: fiducial, sum")->GetXaxis()->SetTitle("First cluster y position (mm)");
    plotter->get2DHistogram("cluster y vs cluster y - cuts: fiducial, sum")->GetYaxis()->SetTitle("Second cluster y position (mm)");

    plotter->build2DHistogram("cluster x vs cluster x - not matched", 100, -200, 100, 100, -200, 100);
    plotter->get2DHistogram("cluster x vs cluster x - not matched")->GetXaxis()->SetTitle("First cluster x position (mm)");
    plotter->get2DHistogram("cluster x vs cluster x - not matched")->GetYaxis()->SetTitle("Second cluster x position (mm)");
    
    plotter->build2DHistogram("cluster y position - not matched", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster y position - not matched")->GetXaxis()->SetTitle("First cluster y position (mm)");
    plotter->get2DHistogram("cluster y position - not matched")->GetYaxis()->SetTitle("Second cluster y position (mm)");

    // Invariant mass
    plotter->build1DHistogram("invariant mass - mollers", 50, 0, 0.1)->GetXaxis()->SetTitle("Invariant mass (GeV)");
    
    plotter->build2DHistogram("p v theta - candidates", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - candidates")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - candidates")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("p v theta - matched", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - matched")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - matched")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("p v theta - not matched", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - not matched")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - not matched")->GetYaxis()->SetTitle("#theta");

}

std::string TagProbeAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;
}

bool TagProbeAnalysis::passClusterEnergySumCut(EcalCluster* first_cluster, EcalCluster* second_cluster) {
    
    double cluster_energy_sum = first_cluster->getEnergy() + second_cluster->getEnergy(); 
    if (cluster_energy_sum > 1.15 || cluster_energy_sum < .75) return false; 

    return true;
}

bool TagProbeAnalysis::passFiducialCut(EcalCluster* first_cluster, EcalCluster* second_cluster) { 
   
    // Require that they are on the electron side in x
    if (first_cluster->getPosition()[0] > 0 || second_cluster->getPosition()[0] > 0) return false;
   
    double cluster_pair_sum_x = first_cluster->getPosition()[0] + second_cluster->getPosition()[0]; 

    if (cluster_pair_sum_x < -180 || cluster_pair_sum_x > -140) return false;

    /* 
    double cluster_pair_delta_x = first_cluster->getPosition()[0] - second_cluster->getPosition()[0]; 

    if (abs(first_cluster->getPosition()[0] - second_cluster->getPosition()[0]) > 60) return false; 
    
    */

    return true;
}
