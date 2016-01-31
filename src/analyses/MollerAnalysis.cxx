/**
 * @file MollerAnalysis.cxx
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date July 10, 2015
 */

#include <MollerAnalysis.h>

MollerAnalysis::MollerAnalysis()
    : plotter(new Plotter()),
      matcher(new TrackClusterMatcher()),
      ecal_utils(new EcalUtils()),
      class_name("MollerAnalysis"), 
      event_counter(0),
      bias_on_counter(0), 
      single1_trigger_counter(0),
      svt_closed_position_counter(0),
      cluster_pair_counter(0) {

}

MollerAnalysis::~MollerAnalysis() {
    delete plotter;
    delete matcher; 
    delete ecal_utils; 
}

void MollerAnalysis::initialize() { 
    this->bookHistograms(); 
    matcher->enablePlots(true);
}

void MollerAnalysis::processEvent(HpsEvent* event) { 
   
    // Increment the total events counter
    event_counter++; 

    // Get a "good" pair from the event.  If a good pair isn't found, skip
    // the event.
    std::vector<EcalCluster*> pair = ecal_utils->getClusterPair(event);
    if (pair[0] == nullptr || pair[1] == nullptr) return;
    cluster_pair_counter++;
 
    // Require that the two clusters are on the electron side
    if (pair[0]->getPosition()[0] > 0 || pair[1]->getPosition()[0] > 0) return; 

    std::vector<double> cluster_energy = { 
        pair[0]->getEnergy(), 
        pair[1]->getEnergy()
    };
    
    double cluster_energy_sum = cluster_energy[0] + cluster_energy[1]; 
    double cluster_energy_diff = cluster_energy[0] - cluster_energy[1];

    std::vector<double> cluster_time = { 
        pair[0]->getClusterTime(), 
        pair[1]->getClusterTime()
    };

    double cluster_time_diff = cluster_time[0] - cluster_time[1]; 

    std::vector<double> cluster_x = {
        pair[0]->getPosition()[0], 
        pair[1]->getPosition()[0] 
    };

    std::vector<double> cluster_y = {
        pair[0]->getPosition()[1], 
        pair[1]->getPosition()[1] 
    };

    double cluster_x_sum = cluster_x[0] + cluster_x[1]; 
    double cluster_x_diff = cluster_x[0] - cluster_x[1]; 

    std::string cuts = "cuts: electron_side";
    
    plotter->get1DHistogram("cluster pair energy sum - " + cuts)->Fill(cluster_energy_sum);
    plotter->get1DHistogram("cluster pair energy diff - " + cuts)->Fill(cluster_energy_diff);
    plotter->get1DHistogram("cluster pair time dt - " + cuts)->Fill(cluster_time_diff);
    plotter->get1DHistogram("cluster pair delta x - " + cuts)->Fill(cluster_x_diff);
    plotter->get1DHistogram("cluster pair x sum - " + cuts)->Fill(cluster_x_sum);
    plotter->get2DHistogram("cluster pair energy - " + cuts)->Fill(cluster_energy[0], cluster_energy[1]);
    plotter->get2DHistogram("cluster pair time - " + cuts)->Fill(cluster_time[0],cluster_time[1]); 
    plotter->get2DHistogram("cluster x vs cluster x - " + cuts)->Fill(cluster_x[0], cluster_x[1]);
    plotter->get2DHistogram("cluster y vs cluster y - " + cuts)->Fill(cluster_y[0], cluster_y[1]);

    // Find all track-cluster matches in the event
    matcher->findAllMatches(event);

    // Check if matches were found for the two clusters. If a match wasn't
    // found, skip the event.
    if (matcher->getMatchingTrack(pair[0]) == nullptr || matcher->getMatchingTrack(pair[1]) == nullptr) return;
    
    std::vector<SvtTrack*> tracks = { 
        matcher->getMatchingTrack(pair[0]), 
        matcher->getMatchingTrack(pair[1])
    };

    // Require that both tracks are negatively charged. If the event otherwise.
    if ((tracks[0]->getCharge() + tracks[1]->getCharge()) != -2) return;

    std::vector<double> track_time = { 
        tracks[0]->getTrackTime(), 
        tracks[1]->getTrackTime() 
    };

    double track_pair_dt = tracks[0]->getTrackTime() - tracks[1]->getTrackTime();     
    
    std::vector<double> p = { 
        AnalysisUtils::getMagnitude(tracks[0]->getMomentum()),
        AnalysisUtils::getMagnitude(tracks[1]->getMomentum()) 
    };
    
    double p_sum = p[0] + p[1];

    std::vector<double> track_theta = { 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(tracks[0]))), 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(tracks[1])))
    };
  
    cuts += ", matched, e-e-";

    plotter->get1DHistogram("cluster pair energy sum - " + cuts)->Fill(cluster_energy_sum);
    plotter->get1DHistogram("cluster pair energy diff - " + cuts)->Fill(cluster_energy_diff);
    plotter->get1DHistogram("cluster pair time dt - " + cuts)->Fill(cluster_time_diff);
    plotter->get1DHistogram("cluster pair delta x - " + cuts)->Fill(cluster_x_diff);
    plotter->get1DHistogram("cluster pair x sum - " + cuts)->Fill(cluster_x_sum);
    plotter->get2DHistogram("cluster pair energy - " + cuts)->Fill(cluster_energy[0], cluster_energy[1]);
    plotter->get2DHistogram("cluster pair time - " + cuts)->Fill(cluster_time[0],cluster_time[1]); 
    plotter->get2DHistogram("cluster x vs cluster x - " + cuts)->Fill(cluster_x[0], cluster_x[1]);
    plotter->get2DHistogram("cluster y vs cluster y - " + cuts)->Fill(cluster_y[0], cluster_y[1]);

    plotter->get1DHistogram("track pair dt - " + cuts)->Fill(track_pair_dt);
    plotter->get1DHistogram("p sum - " + cuts)->Fill(p_sum);
    plotter->get2DHistogram("p[e-] v p[e-] - " + cuts)->Fill(p[0], p[1]);
    plotter->get2DHistogram("track time - " + cuts)->Fill(track_time[0], track_time[1]);
    plotter->get2DHistogram("track theta - " + cuts)->Fill(track_theta[0], track_theta[1]); 
    plotter->get2DHistogram("p v theta - " + cuts)->Fill(p[0], track_theta[0]); 
    plotter->get2DHistogram("p v theta - " + cuts)->Fill(p[1], track_theta[1]); 

    for (auto& track : tracks) { 
        double track_p = AnalysisUtils::getMagnitude(track->getMomentum()); 
        track->isTopTrack() ? plotter->get1DHistogram("p top - " + cuts)->Fill(track_p) 
            : plotter->get1DHistogram("p bottom - " + cuts)->Fill(track_p);
    }

    if (std::abs(track_pair_dt) > 4) return;

    cuts += ", track dt";
    plotter->get1DHistogram("track pair dt - " + cuts)->Fill(track_pair_dt);
    plotter->get1DHistogram("p sum - " + cuts)->Fill(p_sum);
    plotter->get2DHistogram("p[e-] v p[e-] - " + cuts)->Fill(p[0], p[1]);
    plotter->get2DHistogram("track time - " + cuts)->Fill(track_time[0], track_time[1]);
    plotter->get2DHistogram("track theta - " + cuts)->Fill(track_theta[0], track_theta[1]); 
    plotter->get2DHistogram("p v theta - " + cuts)->Fill(p[0], track_theta[0]); 
    plotter->get2DHistogram("p v theta - " + cuts)->Fill(p[1], track_theta[1]); 

    if (p_sum > 1.2) return;
    if (p_sum < .88) return;

    cuts = "moller";

    plotter->get1DHistogram("cluster pair energy sum - " + cuts)->Fill(cluster_energy_sum);
    plotter->get1DHistogram("cluster pair energy diff - " + cuts)->Fill(cluster_energy_diff);
    plotter->get1DHistogram("cluster pair time dt - " + cuts)->Fill(cluster_time_diff);
    plotter->get1DHistogram("cluster pair delta x - " + cuts)->Fill(cluster_x_diff);
    plotter->get1DHistogram("cluster pair x sum - " + cuts)->Fill(cluster_x_sum);
    plotter->get2DHistogram("cluster pair energy - " + cuts)->Fill(cluster_energy[0], cluster_energy[1]);
    plotter->get2DHistogram("cluster pair time - " + cuts)->Fill(cluster_time[0],cluster_time[1]); 
    plotter->get2DHistogram("cluster x vs cluster x - " + cuts)->Fill(cluster_x[0], cluster_x[1]);
    plotter->get2DHistogram("cluster y vs cluster y - " + cuts)->Fill(cluster_y[0], cluster_y[1]);
    plotter->get2DHistogram("cluster position - " + cuts)->Fill(cluster_x[0], cluster_y[0]);
    plotter->get2DHistogram("cluster position - " + cuts)->Fill(cluster_x[1], cluster_y[1]);
    
    plotter->get1DHistogram("track pair dt - " + cuts)->Fill(track_pair_dt);
    plotter->get1DHistogram("p sum - " + cuts)->Fill(p_sum);
    plotter->get2DHistogram("p[e-] v p[e-] - " + cuts)->Fill(p[0], p[1]);
    plotter->get2DHistogram("track time - " + cuts)->Fill(track_time[0], track_time[1]);
    plotter->get2DHistogram("track theta - " + cuts)->Fill(track_theta[0], track_theta[1]); 
    plotter->get2DHistogram("p v theta - " + cuts)->Fill(p[0], track_theta[0]); 
    plotter->get2DHistogram("p v theta - " + cuts)->Fill(p[1], track_theta[1]); 

    double mass = AnalysisUtils::getInvariantMass(tracks[0], tracks[1]); 
    plotter->get1DHistogram("invariant mass - mollers")->Fill(mass);

    std::vector<GblTrack*> gbl_tracks = { 
        nullptr,
        nullptr
    };

    //std::cout << "Total GBL tracks: " << event->getNumberOfGblTracks() << std::endl;
    for (int gbl_track_n = 0; gbl_track_n < event->getNumberOfGblTracks(); ++gbl_track_n) { 
        GblTrack* gbl_track = event->getGblTrack(gbl_track_n); 
        
        if (gbl_track->getSeedTrack() == tracks[0]) { 
            //std::cout << "First GBL track match found." << std::endl;
            gbl_tracks[0] = gbl_track;
        } else if (gbl_track->getSeedTrack() == tracks[1]) { 
            //std::cout << "Second GBL track match found." << std::endl;
            gbl_tracks[1] = gbl_track;
        }
    }

    if (event->getNumberOfParticles(HpsParticle::TC_MOLLER_CANDIDATE) == 0) {
        std::cout << "Moller found but event doesn't contain TC moller candidate." << std::endl;
        return;
    }
    
    for (int particle_n = 0; particle_n < event->getNumberOfParticles(HpsParticle::TC_MOLLER_CANDIDATE); ++particle_n) {
        HpsParticle* particle = event->getParticle(HpsParticle::TC_MOLLER_CANDIDATE, particle_n);

        TRefArray* daughter_particles = particle->getParticles();
        
        if (((HpsParticle*) daughter_particles->At(0))->getType() > 32 && event->getNumberOfParticles(HpsParticle::TC_MOLLER_CANDIDATE) == 2) {
            //std::cout << "Total TC candidates: " << event->getNumberOfParticles(HpsParticle::TC_MOLLER_CANDIDATE) << std::endl;
            if (gbl_tracks[0] == nullptr || gbl_tracks[1] == nullptr) continue;
            //std::cout << "GBL particle found" << std::endl; 
            plotter->get1DHistogram("invariant mass - mollers - GBL TC")->Fill(particle->getMass());
        
            plotter->get1DHistogram("v_{x} - mollers - GBL TC")->Fill(particle->getVertexPosition()[0]); 
            plotter->get1DHistogram("v_{y} - mollers - GBL TC")->Fill(particle->getVertexPosition()[1]); 
            plotter->get1DHistogram("v_{z} - mollers - GBL TC")->Fill(particle->getVertexPosition()[2]); 

            plotter->get2DHistogram("invariant mass : cluster pair x sum - mollers - GBL TC")->Fill(particle->getMass(), cluster_x_sum); 
            plotter->get2DHistogram("invariant mass : cluster pair delta - mollers - GBL TC")->Fill(particle->getMass(), cluster_x_diff);

            plotter->get2DHistogram("invariant mass : Vextex x - mollers - GBL TC")->Fill(particle->getMass(), particle->getVertexPosition()[0]);
            plotter->get2DHistogram("invariant mass : Vextex y - mollers - GBL TC")->Fill(particle->getMass(), particle->getVertexPosition()[1]);
            plotter->get2DHistogram("invariant mass : Vextex #Chi^{2} - mollers - GBL TC")->Fill(particle->getMass(), particle->getVertexFitChi2());

            continue;
        } 
        
        if (daughter_particles->IndexOf(tracks[0]->getParticle())
                *daughter_particles->IndexOf(tracks[1]->getParticle()) != 0) continue;        

        //std::cout << "Seed particle found" << std::endl;
        plotter->get1DHistogram("invariant mass - mollers - TC")->Fill(particle->getMass());
    
        plotter->get1DHistogram("v_{x} - mollers - TC")->Fill(particle->getVertexPosition()[0]); 
        plotter->get1DHistogram("v_{y} - mollers - TC")->Fill(particle->getVertexPosition()[1]); 
        plotter->get1DHistogram("v_{z} - mollers - TC")->Fill(particle->getVertexPosition()[2]); 

        plotter->get2DHistogram("invariant mass : cluster pair x sum - mollers - TC")->Fill(particle->getMass(), cluster_x_sum); 
        plotter->get2DHistogram("invariant mass : cluster pair delta - mollers - TC")->Fill(particle->getMass(), cluster_x_diff);

        plotter->get2DHistogram("invariant mass : Vextex x - mollers - TC")->Fill(particle->getMass(), particle->getVertexPosition()[0]);
        plotter->get2DHistogram("invariant mass : Vextex y - mollers - TC")->Fill(particle->getMass(), particle->getVertexPosition()[1]);
        plotter->get2DHistogram("invariant mass : Vextex #Chi^{2} - mollers - TC")->Fill(particle->getMass(), particle->getVertexFitChi2());
    
    }

    /*

    if (gbl_tracks[0] == nullptr || gbl_tracks[1] == nullptr) return;

    */
}

void MollerAnalysis::finalize() {

    ecal_utils->saveHistograms();
    matcher->saveHistograms();

    plotter->saveToPdf("moller_analysis.pdf");
    plotter->saveToRootFile("moller_analysis.root");
}


void MollerAnalysis::bookHistograms() {

    TH1* plot = nullptr; 

    std::vector<std::string> cuts = { 
        "cuts: electron_side", 
        "cuts: electron_side, matched, e-e-",
        "cuts: electron_side, matched, e-e-, track dt",
        "moller"
    };
    
    for (auto& cut : cuts) { 
        
        //----------//
        //   Ecal   //
        //----------//
         
        plot = plotter->build1DHistogram("cluster pair energy sum - " + cut, 100, 0, 1.5);
        plot->GetXaxis()->SetTitle("Cluster pair energy sum (GeV)");

        plot = plotter->build1DHistogram("cluster pair energy diff - " + cut, 100, -1, 1);
        plot->GetXaxis()->SetTitle("Cluster pair energy difference (GeV)");
    
        plot = plotter->build1DHistogram("cluster pair time dt - " + cut, 40, -10, 10);
    
        plot = plotter->build1DHistogram("cluster pair delta x - " + cut, 100, -200, 200);
        plot->GetXaxis()->SetTitle("Ecal cluster pair dx (mm)");
    
        plot = plotter->build1DHistogram("cluster pair x sum - " + cut, 100, -250, -100);
        plot->GetXaxis()->SetTitle("Ecal cluster pair x sum (mm)");
        
        plot = plotter->build2DHistogram("cluster pair energy - " + cut, 50, 0, 1.5, 50, 0, 1.5);
        plot->GetXaxis()->SetTitle("Cluster energy (GeV)");
        plot->GetYaxis()->SetTitle("Cluster energy (GeV)");
        
        plot = plotter->build2DHistogram("cluster pair time - " + cut, 250, 0, 125, 250, 0, 125);
        plot->GetXaxis()->SetTitle("Cluster time (ns)");
        plot->GetYaxis()->SetTitle("Cluster time (ns)");
    
        plot = plotter->build2DHistogram("cluster x vs cluster x - " + cut, 200, -200, 200, 200, -200, 200);
        plot->GetXaxis()->SetTitle("Cluster position - x (mm)");
        plot->GetYaxis()->SetTitle("Cluster position - x (mm)");
    
        plot = plotter->build2DHistogram("cluster y vs cluster y - " + cut, 200, -200, 200, 200, -200, 200);
        plot->GetYaxis()->SetTitle("Cluster position - y (mm)");
        plot->GetYaxis()->SetTitle("Cluster position - y (mm)"); 

        if (cut.compare(cuts[0]) == 0) continue; 
    
        plot = plotter->build1DHistogram("track pair dt - " + cut, 100, -10, 10);
        plot->GetXaxis()->SetTitle("#Delta Track time [ns]");
        
        plotter->build1DHistogram("p sum - " + cut, 100, 0, 2.0)->GetXaxis()->SetTitle("p sum (GeV)");
   
        plot = plotter->build2DHistogram("p[e-] v p[e-] - " + cut, 50, 0, 1.5, 50, 0, 1.5);
        plot->GetXaxis()->SetTitle("p[e-] [GeV]"); 
        plot->GetYaxis()->SetTitle("p[e-] [GeV]"); 

        plot = plotter->build2DHistogram("track time - " + cut, 100, -10, 10, 100, -10, 10);
        plot->GetXaxis()->SetTitle("Track time [ns]");
        plot->GetYaxis()->SetTitle("Track time [ns]");

        plot = plotter->build2DHistogram("track theta - " + cut, 40, 0, 0.1, 40, 0, 0.1);
        plot->GetXaxis()->SetTitle("#theta");
        plot->GetYaxis()->SetTitle("#theta");

        plot = plotter->build2DHistogram("p v theta - " + cut, 50, 0, 1.5, 40, 0, 0.1);
        plot->GetXaxis()->SetTitle("p (GeV)");
        plot->GetYaxis()->SetTitle("#theta");

        plotter->build1DHistogram("p top - " + cut, 50, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");
        plotter->build1DHistogram("p bottom - " + cut, 50, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");

    }
   

    // Invariant mass
    plot = plotter->build1DHistogram("invariant mass - mollers", 100, 0.02, 0.05);
    plot->GetXaxis()->SetTitle("Invariant mass (GeV)");

    plot = plotter->build1DHistogram("invariant mass - mollers - TC", 100, 0.02, 0.05); 
    plot->GetXaxis()->SetTitle("Invariant mass (GeV)"); 

    plot = plotter->build1DHistogram("invariant mass - mollers - GBL TC", 100, 0.02, 0.05);
    plot->GetXaxis()->SetTitle("Invariant mass (GeV)");

    plot = plotter->build1DHistogram("v_{x} - mollers - TC", 100, -1, 1); 
    plot->GetXaxis()->SetTitle("Vertex x position (mm)"); 

    plot = plotter->build1DHistogram("v_{y} - mollers - TC", 100, -1, 1); 
    plot->GetXaxis()->SetTitle("Vertex y position (mm)"); 

    plot = plotter->build1DHistogram("v_{z} - mollers - TC", 100, -1, 1); 
    plot->GetXaxis()->SetTitle("Vertex z position (mm)"); 

    plot = plotter->build1DHistogram("v_{z} - mollers - GBL TC", 100, -1, 1); 
    plot->GetXaxis()->SetTitle("Vertex z position (mm)"); 

    plot = plotter->build1DHistogram("v_{x} - mollers - GBL TC", 100, -1, 1); 
    plot->GetXaxis()->SetTitle("Vertex x position (mm)"); 

    plot = plotter->build1DHistogram("v_{y} - mollers - GBL TC", 100, -1, 1); 
    plot->GetXaxis()->SetTitle("Vertex y position (mm)"); 

    plotter->build2DHistogram("invariant mass : cluster pair x sum - mollers - TC", 100, 0, 0.1, 100, -250, -100);
    plotter->build2DHistogram("invariant mass : cluster pair delta - mollers - TC", 100, 0, 0.1, 100, -150, 200);

    plotter->build2DHistogram("invariant mass : Vextex x - mollers - TC", 100, 0, 0.1, 100, -1, 1); 

    plotter->build2DHistogram("invariant mass : Vextex y - mollers - TC", 100, 0, 0.1, 100, -1, 1); 
    
    plotter->build2DHistogram("invariant mass : Vextex #Chi^{2} - mollers - TC", 100, 0, 0.1, 100, 0, 50);

    plotter->build2DHistogram("invariant mass : cluster pair x sum - mollers - GBL TC", 100, 0, 0.1, 100, -250, -100);
    plotter->build2DHistogram("invariant mass : cluster pair delta - mollers - GBL TC", 100, 0, 0.1, 100, -150, 200);

    plotter->build2DHistogram("invariant mass : Vextex x - mollers - GBL TC", 100, 0, 0.1, 100, -1, 1); 

    plotter->build2DHistogram("invariant mass : Vextex y - mollers - GBL TC", 100, 0, 0.1, 100, -1, 1); 
    
    plotter->build2DHistogram("invariant mass : Vextex #Chi^{2} - mollers - GBL TC", 100, 0, 0.1, 100, 0, 50);

    plotter->build2DHistogram("cluster position - moller", 200, -200, 200, 100, -100, 100);
} 

std::string MollerAnalysis::toString() { 

    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
