/**
 *
 * @file TridentAnalysis.cxx
 * @brief Analysis used to look at Tridents.
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date September 29, 2015
 *
 */

#include <TridentAnalysis.h>

TridentAnalysis::TridentAnalysis()
    : plotter(new Plotter()), 
      ecal_utils(new EcalUtils()), 
      matcher(new TrackClusterMatcher()), 
      class_name("TridentAnalysis"), 
      event_counter(0), 
      good_cluster_pair_event_counter(0),  
      matched_event_counter(0), 
      v0_cand_counter(0) { 
}

TridentAnalysis::~TridentAnalysis() { 
    delete plotter;
    delete ecal_utils; 
    delete matcher;
}

void TridentAnalysis::initialize() { 
    this->bookHistograms(); 
}

void TridentAnalysis::processEvent(HpsEvent* event) { 
   
    event_counter++; 

    // Get a "good" pair from the event.  If a good pair isn't found, skip
    // the event.
    std::vector<EcalCluster*> pair = ecal_utils->getClusterPair(event);
    if (pair.size() != 2 || pair[0] == nullptr || pair[1] == nullptr) return;
    good_cluster_pair_event_counter++;

    // Find all track-cluster matches in the event
    matcher->findAllMatches(event);

    // Check if matches for the two clusters were found
    if (matcher->getMatchingTrack(pair[0]) == nullptr || matcher->getMatchingTrack(pair[1]) == nullptr) return;
    matched_event_counter++; 

    //
    // Define some commonly used variables
    //
   
    // Sort the clusters in order of highest energy 
    int first_cluster_index = 0; int second_cluster_index = 1;
    if (pair[0]->getEnergy() < pair[1]->getEnergy()) { first_cluster_index = 1; second_cluster_index = 0; } 
    std::vector<double> cluster_energy = { 
        pair[first_cluster_index]->getEnergy(), 
        pair[second_cluster_index]->getEnergy()
    };
    
    double cluster_energy_sum = cluster_energy[0] + cluster_energy[1]; 
    double cluster_energy_diff = cluster_energy[0] - cluster_energy[1];

    std::vector<double> cluster_time = { 
        pair[first_cluster_index]->getClusterTime(), 
        pair[second_cluster_index]->getClusterTime()
    };

    //
    // Fill plots of tracks and clusters that are matched to each other
    //

    std::vector<SvtTrack*> tracks = { 
        matcher->getMatchingTrack(pair[first_cluster_index]), 
        matcher->getMatchingTrack(pair[second_cluster_index])
    };

    plotter->get2DHistogram("cluster pair energy - matched")->Fill(cluster_energy[0], cluster_energy[1]); 
    plotter->get1DHistogram("cluster pair energy sum - matched")->Fill(cluster_energy_sum); 
    plotter->get1DHistogram("cluster pair energy diff - matched")->Fill(cluster_energy_diff);
    plotter->get1DHistogram("track time - cluster time - matched")->Fill(tracks[0]->getTrackTime() - cluster_time[0]);  
    plotter->get1DHistogram("track time - cluster time - matched")->Fill(tracks[1]->getTrackTime() - cluster_time[1]);  
    plotter->get1DHistogram("track pair dt - matched")->Fill(tracks[0]->getTrackTime() - tracks[1]->getTrackTime()); 

    // Skip events where both tracks have the same charge
    if (tracks[0]->getCharge()*tracks[1]->getCharge() > 0) return; 

    // Identify the electron and positron in the event
    SvtTrack* electron = tracks[0];
    SvtTrack* positron = tracks[1];
    if (tracks[1]->getCharge() == -1) { electron = tracks[1]; positron = tracks[0]; } 
    v0_cand_counter++;

    // Calculate the momentum of the electron and positrons 
    double electron_p = AnalysisUtils::getMagnitude(electron->getMomentum());
    double positron_p = AnalysisUtils::getMagnitude(positron->getMomentum());

    //
    // Fill plots of events that have e+e- tracks matched to clusters
    //
    plotter->get2DHistogram("cluster pair energy - matched, e+e-")->Fill(cluster_energy[0], cluster_energy[1]); 
    plotter->get1DHistogram("cluster pair energy sum - matched, e+e-")->Fill(cluster_energy_sum); 
    plotter->get1DHistogram("cluster pair energy diff - matched, e+e-")->Fill(cluster_energy_diff);
    plotter->get1DHistogram("track time - cluster time - matched, e+e-")->Fill(tracks[0]->getTrackTime()
                                                                                - cluster_time[0]);  
    plotter->get1DHistogram("track time - cluster time - matched, e+e-")->Fill(tracks[1]->getTrackTime() 
                                                                                - cluster_time[1]);  
    plotter->get1DHistogram("track pair dt - matched, e+e-")->Fill(tracks[0]->getTrackTime() 
                                                                    - tracks[1]->getTrackTime()); 
    plotter->get2DHistogram("p[e+] v p[e-] - matched, e+e-")->Fill(electron_p, positron_p);  
    plotter->get1DHistogram("invariant mass")->Fill(AnalysisUtils::getInvariantMass(electron, positron)); 

    //
    // Get the reconstructed particle associated with both tracks
    //

    if (event->getNumberOfParticles(HpsParticle::UC_V0_CANDIDATE) == 0) {
        std::cout << "Event doesn't contain uncontrained candidate but the event contains an e+e-!" << std::endl;
        return;
    }

    for (int particle_n = 0; particle_n < event->getNumberOfParticles(HpsParticle::UC_V0_CANDIDATE); ++particle_n) {
        
        HpsParticle* particle = event->getParticle(HpsParticle::UC_V0_CANDIDATE, particle_n);
        
        TRefArray* daughter_particles = particle->getParticles();

        if (daughter_particles->FindObject(electron) == 0 || daughter_particles->FindObject(positron) == 0) continue;
    
        std::cout << "Matching candidate has been found" << std::endl;
    }  
}

void TridentAnalysis::finalize() { 
    ecal_utils->saveHistograms();
    plotter->saveToRootFile("trident_analysis.root");

    std::cout << "Trident's passing cuts: " << v0_cand_counter << "/" << event_counter << " = " 
              << v0_cand_counter/event_counter << " %" << std::endl; 
}

void TridentAnalysis::bookHistograms() { 

    // Matched Ecal clusters
    plotter->build1DHistogram("cluster pair energy diff - matched", 150, -1, 1)->GetXaxis()->SetTitle(
            "Cluster Pair Energy Difference (GeV)");
    
    plotter->build1DHistogram("track time - cluster time - matched", 60, 30, 60)->GetXaxis()->SetTitle(
            "Track time - Ecal cluster time"); 

    plotter->build2DHistogram("cluster pair energy - matched", 150, 0, 1.5, 150, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - matched")->GetXaxis()->SetTitle("Cluster Energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - matched")->GetYaxis()->SetTitle("Cluster Energy (GeV)");

    plotter->build1DHistogram("cluster pair energy sum - matched", 150, 0, 1.5)->GetXaxis()->SetTitle(
            "Cluster Pair Energy Sum (GeV)");
    
    plotter->build1DHistogram("cluster pair energy diff - matched, e+e-", 150, -1, 1)->GetXaxis()->SetTitle(
            "Cluster Pair Energy Difference (GeV)");

    plotter->build1DHistogram("track time - cluster time - matched, e+e-", 60, 30, 60)->GetXaxis()->SetTitle(
            "Track time - Ecal cluster time");

    plotter->build1DHistogram("cluster pair energy sum - matched, e+e-", 150, 0, 1.5)->GetXaxis()->SetTitle(
            "Cluster Pair Energy Sum (GeV)");

    plotter->build2DHistogram("cluster pair energy - matched, e+e-", 150, 0, 1.5, 150, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - matched, e+e-")->GetXaxis()->SetTitle("Cluster Energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - matched, e+e-")->GetYaxis()->SetTitle("Cluster Energy (GeV)");

    // Matched tracks
    plotter->build1DHistogram("track pair dt - matched", 50, -25, 25)->GetXaxis()->SetTitle("Track time #Delta t");  
    plotter->build1DHistogram("track pair dt - matched, e+e-", 50, -25, 25)->GetXaxis()->SetTitle(
            "Track time #Delta t");  
    plotter->build1DHistogram("p top - matched, e+e-", 150, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");
    plotter->build1DHistogram("p bottom - matched, e+e-", 150, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");

    plotter->build2DHistogram("p[e+] v p[e-] - matched, e+e-", 150, 0, 1.5, 150, 0, 1.5);
    plotter->get2DHistogram("p[e+] v p[e-] - matched, e+e-")->GetXaxis()->SetTitle("p[e-] (GeV)"); 
    plotter->get2DHistogram("p[e+] v p[e-] - matched, e+e-")->GetYaxis()->SetTitle("p[e+] (GeV)"); 

    // Invariant mass
    plotter->build1DHistogram("invariant mass", 200, 0, 0.1)->GetXaxis()->SetTitle("Mass (GeV)");

}

std::string TridentAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
