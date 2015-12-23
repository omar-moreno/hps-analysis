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
      good_cluster_pair_counter(0),  
      matched_event_counter(0), 
      v0_cand_counter(0) { 
}

TridentAnalysis::~TridentAnalysis() { 
    delete plotter;
    delete ecal_utils; 
    delete matcher;
}

void TridentAnalysis::initialize() { 

    // Enable track-cluster matching plots
    matcher->enablePlots(); 

    this->bookHistograms(); 
}

void TridentAnalysis::processEvent(HpsEvent* event) { 
   
    event_counter++; 

    // Get a "good" pair from the event.  If a good pair isn't found, skip
    // the event.
    std::vector<EcalCluster*> pair = ecal_utils->getClusterPair(event);
    if (pair[0] == nullptr || pair[1] == nullptr) return;
    good_cluster_pair_counter++;

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

    std::vector<double> cluster_time = { 
        pair[first_cluster_index]->getClusterTime(), 
        pair[second_cluster_index]->getClusterTime()
    };

    std::vector<double> cluster_x = {
        pair[first_cluster_index]->getPosition()[0], 
        pair[second_cluster_index]->getPosition()[0] 
    };

    std::vector<double> cluster_y = {
        pair[first_cluster_index]->getPosition()[1], 
        pair[second_cluster_index]->getPosition()[1] 
    };

    double cluster_energy_sum = cluster_energy[0] + cluster_energy[1]; 
    double cluster_energy_diff = cluster_energy[0] - cluster_energy[1];

    //
    // Fill histograms
    //

    plotter->get1DHistogram("cluster energy high")->Fill(cluster_energy[0]); 
    plotter->get1DHistogram("cluster energy low")->Fill(cluster_energy[1]); 
    plotter->get1DHistogram("cluster time")->Fill(cluster_time[0]); 
    plotter->get1DHistogram("cluster time")->Fill(cluster_time[1]); 
    plotter->get2DHistogram("cluster position")->Fill(cluster_x[0], cluster_y[0]); 
    plotter->get2DHistogram("cluster position")->Fill(cluster_x[1], cluster_y[1]); 

    // Find all track-cluster matches in the event
    matcher->findAllMatches(event);

    // Check if matches for the two clusters were found
    if (matcher->getMatchingTrack(pair[0]) == nullptr || matcher->getMatchingTrack(pair[1]) == nullptr) return;
    matched_event_counter++; 

    std::vector<SvtTrack*> tracks = { 
        matcher->getMatchingTrack(pair[first_cluster_index]), 
        matcher->getMatchingTrack(pair[second_cluster_index])
    };

    std::vector<double> track_p = { 
        AnalysisUtils::getMagnitude(tracks[0]->getMomentum()),
        AnalysisUtils::getMagnitude(tracks[1]->getMomentum())
    };

    //
    // Fill plots of tracks and clusters that are matched to each other
    //

    plotter->get1DHistogram("cluster energy high - matched")->Fill(cluster_energy[0]); 
    plotter->get1DHistogram("cluster energy low - matched")->Fill(cluster_energy[1]); 
    plotter->get1DHistogram("cluster time - matched")->Fill(cluster_time[0]); 
    plotter->get1DHistogram("cluster time - matched")->Fill(cluster_time[1]); 
    plotter->get1DHistogram("cluster pair energy sum - matched")->Fill(cluster_energy_sum); 
    plotter->get1DHistogram("cluster pair energy diff - matched")->Fill(cluster_energy_diff);
    plotter->get1DHistogram("track time - cluster time - matched")->Fill(tracks[0]->getTrackTime() - cluster_time[0]);  
    plotter->get1DHistogram("track time - cluster time - matched")->Fill(tracks[1]->getTrackTime() - cluster_time[1]);  
    plotter->get1DHistogram("track pair dt - matched")->Fill(tracks[0]->getTrackTime() - tracks[1]->getTrackTime()); 
    plotter->get1DHistogram("p high - matched")->Fill(track_p[0]); 
    plotter->get1DHistogram("p low - matched")->Fill(track_p[1]); 
    plotter->get1DHistogram("p sum - matched")->Fill(track_p[0] + track_p[1]); 
    plotter->get2DHistogram("cluster pair energy - matched")->Fill(cluster_energy[0], cluster_energy[1]); 
    plotter->get2DHistogram("cluster position - matched")->Fill(cluster_x[0], cluster_y[0]); 
    plotter->get2DHistogram("cluster position - matched")->Fill(cluster_x[1], cluster_y[1]); 
    plotter->get2DHistogram("p_{e} v p_{e} - matched")->Fill(track_p[0], track_p[1]);  

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
    plotter->get1DHistogram("cluster pair energy sum - matched, e+e-")->Fill(cluster_energy_sum); 
    plotter->get1DHistogram("cluster pair energy diff - matched, e+e-")->Fill(cluster_energy_diff);
    plotter->get1DHistogram("track time - cluster time - matched, e+e-")->Fill(
            tracks[0]->getTrackTime() - cluster_time[0]);  
    plotter->get1DHistogram("track time - cluster time - matched, e+e-")->Fill(
            tracks[1]->getTrackTime() - cluster_time[1]);  
    plotter->get1DHistogram("track pair dt - matched, e+e-")->Fill(
            tracks[0]->getTrackTime() - tracks[1]->getTrackTime()); 
    plotter->get1DHistogram("invariant mass")->Fill(AnalysisUtils::getInvariantMass(electron, positron)); 
    plotter->get1DHistogram("p high - matched, e+e-")->Fill(track_p[0]); 
    plotter->get1DHistogram("p low - matched, e+e-")->Fill(track_p[1]); 
    plotter->get1DHistogram("p sum - matched, e+e-")->Fill(track_p[0] + track_p[1]); 
    plotter->get2DHistogram("cluster pair energy - matched, e+e-")->Fill(cluster_energy[0], cluster_energy[1]); 
    plotter->get2DHistogram("cluster position - matched, e+e-")->Fill(cluster_x[0], cluster_y[0]); 
    plotter->get2DHistogram("cluster position - matched, e+e-")->Fill(cluster_x[1], cluster_y[1]); 
    plotter->get2DHistogram("p_{e+} v p_{e-} - matched, e+e-")->Fill(electron_p, positron_p);  
    plotter->get2DHistogram("invariant mass v track p sum")->Fill(
            AnalysisUtils::getInvariantMass(electron, positron), track_p[0]+track_p[1]);
    plotter->get2DHistogram("cluster x high v track p sum - matched, e+e-")->Fill(cluster_x[0], track_p[0]+track_p[1]);
    plotter->get2DHistogram("cluster y high v track p sum - matched, e+e-")->Fill(cluster_y[0], track_p[0]+track_p[1]);
    plotter->get2DHistogram("cluster x low v track p sum - matched, e+e-")->Fill(cluster_x[1], track_p[0]+track_p[1]);
    plotter->get2DHistogram("cluster y low v track p sum - matched, e+e-")->Fill(cluster_y[1], track_p[0]+track_p[1]);


    if (track_p[0]+track_p[1] >= 0.8) { 
        plotter->get1DHistogram("invariant mass - final")->Fill(AnalysisUtils::getInvariantMass(electron, positron)); 
    }

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
    
    // Save the track-Ecal cluster matching plots to a ROOT file
    matcher->saveHistograms();  
    ecal_utils->saveHistograms();
    plotter->saveToRootFile("trident_analysis.root");

    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    std::cout << "Number of pair1 triggers: " << event_counter << std::endl;
    std::cout << "Good clusters: " << good_cluster_pair_counter << "/" <<  event_counter << " = "
              << good_cluster_pair_counter/event_counter << " %" << std::endl;
    std::cout << "Trident's passing cuts: " << v0_cand_counter << "/" << event_counter << " = " 
              << v0_cand_counter/event_counter << " %" << std::endl; 
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
}

void TridentAnalysis::bookHistograms() { 

    TH1* plot = nullptr; 
    
    //
    // Plots of good clusters
    //
    plotter->build1DHistogram("cluster energy high", 150, 0, 1.5)->GetXaxis()->SetTitle("Ecal cluster energy (GeV)");
    plotter->build1DHistogram("cluster energy low", 150, 0, 1.5)->GetXaxis()->SetTitle("Ecal cluster energy (GeV)");

    plotter->build1DHistogram("cluster time", 200, 0, 100)->GetXaxis()->SetTitle("Ecal cluster time (ns)");

    plot = plotter->build2DHistogram("cluster position", 200, -200, 200, 100, -100, 100);
    plot->GetXaxis()->SetTitle("Ecal cluster x (mm)");
    plot->GetXaxis()->SetTitle("Ecal cluster y (mm)");
    
    //
    // Plots of matched tracks and clusters   
    //
   
    plot = plotter->build1DHistogram("cluster energy high - matched", 150, 0, 1.5);
    plot->GetXaxis()->SetTitle("Ecal Cluster Energy (GeV)"); 
    
    plot = plotter->build1DHistogram("cluster energy low - matched", 150, 0, 1.5);
    plot->GetXaxis()->SetTitle("Ecal Cluster Energy (GeV)"); 

    plotter->build1DHistogram("cluster time - matched", 200, 0, 100)->GetXaxis()->SetTitle("Ecal cluster time (ns)");

    plot = plotter->build1DHistogram("cluster pair energy diff - matched", 150, -1, 1); 
    plot->GetXaxis()->SetTitle("Cluster Pair Energy Difference (GeV)");

    plot = plotter->build1DHistogram("cluster pair energy sum - matched", 150, 0, 1.5);
    plot->GetXaxis()->SetTitle("Cluster Pair Energy Sum (GeV)");
    
    plot = plotter->build2DHistogram("cluster pair energy - matched", 150, 0, 1.5, 150, 0, 1.5);
    plot->GetXaxis()->SetTitle("Cluster Energy (GeV)");
    plot->GetYaxis()->SetTitle("Cluster Energy (GeV)");

    plot = plotter->build2DHistogram("cluster position - matched", 200, -200, 200, 100, -100, 100);
    plot->GetXaxis()->SetTitle("Ecal cluster x (mm)");
    plot->GetXaxis()->SetTitle("Ecal cluster y (mm)");

    plot = plotter->build1DHistogram("track time - cluster time - matched", 60, 30, 60); 
    plot->GetXaxis()->SetTitle("Track time - Ecal cluster time"); 

    plot = plotter->build1DHistogram("track pair dt - matched", 50, -25, 25);
    plot->GetXaxis()->SetTitle("Track time #Delta t");  

    plotter->build1DHistogram("p high - matched", 150, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");
    plotter->build1DHistogram("p low - matched", 150, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");
    plotter->build1DHistogram("p sum - matched", 150, 0, 1.5)->GetXaxis()->SetTitle("p sum (GeV)"); 

    plot = plotter->build2DHistogram("p_{e} v p_{e} - matched", 150, 0, 1.5, 150, 0, 1.5);
    plot->GetXaxis()->SetTitle("p_{e} (GeV)"); 
    plot->GetYaxis()->SetTitle("p_{e} (GeV)"); 

    //
    // Plots of tracks and clusters where an e+e- has been identified
    //

    plot = plotter->build1DHistogram("cluster pair energy diff - matched, e+e-", 150, -1, 1); 
    plot->GetXaxis()->SetTitle("Cluster Pair Energy Difference (GeV)");

    plotter->build1DHistogram("cluster pair energy sum - matched, e+e-", 150, 0, 1.5); 
    plot->GetXaxis()->SetTitle("Cluster Pair Energy Sum (GeV)");

    plot = plotter->build2DHistogram("cluster pair energy - matched, e+e-", 150, 0, 1.5, 150, 0, 1.5);
    plot->GetXaxis()->SetTitle("Cluster Energy (GeV)");
    plot->GetYaxis()->SetTitle("Cluster Energy (GeV)");

    plot = plotter->build2DHistogram("cluster position - matched, e+e-", 200, -200, 200, 100, -100, 100);
    plot->GetXaxis()->SetTitle("Ecal cluster x (mm)");
    plot->GetXaxis()->SetTitle("Ecal cluster y (mm)");

    plot = plotter->build1DHistogram("track time - cluster time - matched, e+e-", 60, 30, 60);
    plot->GetXaxis()->SetTitle("Track time - Ecal cluster time");

    plot = plotter->build1DHistogram("track pair dt - matched, e+e-", 50, -25, 25);
    plot->GetXaxis()->SetTitle("Track time #Delta t");  


    plotter->build1DHistogram("p high - matched, e+e-", 150, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");
    plotter->build1DHistogram("p low - matched, e+e-", 150, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");
    plotter->build1DHistogram("p sum - matched, e+e-", 150, 0, 1.5)->GetXaxis()->SetTitle("p sum (GeV)"); 
    plotter->build1DHistogram("p sum - matched, e+e-, fiducial", 150, 0, 1.5)->GetXaxis()->SetTitle("p sum (GeV)"); 
    

    plot = plotter->build2DHistogram("p_{e+} v p_{e-} - matched, e+e-", 150, 0, 1.5, 150, 0, 1.5);
    plot->GetXaxis()->SetTitle("p_{e-} (GeV)"); 
    plot->GetYaxis()->SetTitle("p_{e+} (GeV)");

    plot = plotter->build2DHistogram("cluster x high v track p sum - matched, e+e-", 200, -200, 200, 150, 0, 1.5);
    plot->GetXaxis()->SetTitle("Ecal cluster x (mm)");
    plot->GetYaxis()->SetTitle("p Sum (GeV)"); 

    plot = plotter->build2DHistogram("cluster x low v track p sum - matched, e+e-", 200, -200, 200, 150, 0, 1.5);
    plot->GetXaxis()->SetTitle("Ecal cluster x (mm)");
    plot->GetYaxis()->SetTitle("p Sum (GeV)"); 

    plot = plotter->build2DHistogram("cluster y high v track p sum - matched, e+e-", 100, -100, 100, 150, 0, 1.5);
    plot->GetXaxis()->SetTitle("Ecal cluster y (mm)");
    plot->GetYaxis()->SetTitle("p Sum (GeV)"); 

    plot = plotter->build2DHistogram("cluster y low v track p sum - matched, e+e-", 100, -100, 100, 150, 0, 1.5);
    plot->GetXaxis()->SetTitle("Ecal cluster y (mm)");
    plot->GetYaxis()->SetTitle("p Sum (GeV)"); 

    //plot = plotter->build1DHistogram("cos(#theta_{e+}) v cos(#theta_{e-}) - matched, e+e-");
    //plot->GetXaxis()->SetTitle("cos(#theta_{e-})");  
    //plot->GetYaxis()->SetTitle("cos(#theta_{e-})");  

    //
    // Invariant mass
    //
    
    plotter->build1DHistogram("invariant mass", 1000, 0, 0.1)->GetXaxis()->SetTitle("Mass (GeV)");
    plotter->build1DHistogram("invariant mass - final", 1000, 0, 0.1)->GetXaxis()->SetTitle("Mass (GeV)");

    plot = plotter->build2DHistogram("invariant mass v track p sum", 200, 0, 0.1, 150, 0, 1.5);
    plot->GetXaxis()->SetTitle("Mass (GeV)");
    plot->GetYaxis()->SetTitle("p Sum (GeV)"); 
}

std::string TridentAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
