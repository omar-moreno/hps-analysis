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

    // Require that both tracks are negatively charged. If the event otherwise.
    if ((tracks[0]->getCharge() + tracks[1]->getCharge()) != -2) return;

    std::vector<SvtTrack*> tracks = { 
        matcher->getMatchingTrack(pair[0]), 
        matcher->getMatchingTrack(pair[1])
    };

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
    plotter->get2DHistogram("p[e-] v p[e-] - matched, time")->Fill(p[0], p[1]);
    plotter->get2DHistogram("track time - " + cuts)->Fill(track_time[0], track_time[1]);
    plotter->get2DHistogram("track theta - " + cuts)->Fill(track_theta[0], track_theta[1]); 
    plotter->get2DHistogram("p v theta - " + cuts)->Fill(p[0], track_theta[0]); 
    plotter->get2DHistogram("p v theta - " + cuts)->Fill(p[1], track_theta[1]); 

    if (p_sum > 1.2) return;
    if (p_sum < .88) return;

    plotter->get1DHistogram("cluster pair energy sum - moller")->Fill(cluster_energy_sum);
    plotter->get1DHistogram("cluster pair energy diff - moller")->Fill(cluster_energy_diff);
    plotter->get1DHistogram("cluster pair delta x - moller")->Fill(cluster_x_diff);
    plotter->get1DHistogram("cluster pair x sum - moller")->Fill(cluster_x_sum);
    plotter->get2DHistogram("cluster pair energy - moller")->Fill(cluster_energy[0], cluster_energy[1]);
    plotter->get2DHistogram("cluster x vs cluster x - moller")->Fill(cluster_x[0], cluster_x[1]);
    plotter->get2DHistogram("cluster y vs cluster y - moller")->Fill(cluster_y[0], cluster_y[1]);
    plotter->get2DHistogram("cluster position - moller")->Fill(cluster_x[0], cluster_y[0]);
    plotter->get2DHistogram("cluster position - moller")->Fill(cluster_x[1], cluster_y[1]);
    plotter->get2DHistogram("p[e-] v p[e-] - moller")->Fill(p[0], p[1]);
    plotter->get2DHistogram("track theta - moller")->Fill(track_theta[0], track_theta[1]); 
    plotter->get2DHistogram("p v theta - moller")->Fill(p[0], track_theta[0]); 
    plotter->get2DHistogram("p v theta - moller")->Fill(p[1], track_theta[1]); 

    double mass = AnalysisUtils::getInvariantMass(tracks[0], tracks[1]); 
    plotter->get1DHistogram("invariant mass - mollers")->Fill(mass);

    std::vector<GblTrack*> gbl_tracks = { 
        nullptr, 
        nullptr
    };

    for (int gbl_track_n = 0; gbl_track_n < event->getNumberOfGblTracks(); ++gbl_track_n) { 
        GblTrack* gbl_track = event->getGblTrack(gbl_track_n); 

        if (gbl_track->getSeedTrack() == tracks[0]) { 
            gbl_tracks[0] = gbl_track;
        } else if (gbl_track->getSeedTrack() == tracks[1]) { 
            gbl_tracks[1] = gbl_track;
        }
    }

    if (gbl_tracks[0] == nullptr || gbl_tracks[1] == nullptr) return;

    mass = AnalysisUtils::getInvariantMass(tracks[0], tracks[1]); 
    plotter->get1DHistogram("invariant mass - mollers - GBL")->Fill(mass);
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
         
        plot = plotter->build1DHistogram("cluster pair energy sum - " + cut, 50, 0, 1.5)
        plot->GetXaxis()->SetTitle("Cluster pair energy sum (GeV)");

        plot = plotter->build1DHistogram("cluster pair energy diff - " + cut, 50, -1, 1);
        plot->GetXaxis()->SetTitle("Cluster pair energy difference (GeV)");
    
        plot = plotter->build1DHistogram("cluster pair time dt - " + cut, 40, -10, 10);
    
        plot = plotter->build1DHistogram("cluster pair delta x - " + cut, 50, -200, 200);
        plot->GetXaxis()->SetTitle("Ecal cluster pair dx (mm)");
    
        plot = plotter->build1DHistogram("cluster pair x sum - " + cut, 50, -250, -100);
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
    }
   
    plotter->build1DHistogram("p sum - matched, e-e-, time", 100, 0, 2.0)->GetXaxis()->SetTitle("p sum (GeV)");
   
    plotter->build1DHistogram("p top - matched", 50, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");
    plotter->build1DHistogram("p bottom - matched", 50, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");

    plotter->build2DHistogram("p[e-] v p[e-] - matched, time", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e-] v p[e-] - matched, time")->GetXaxis()->SetTitle("p[e-] [GeV]"); 
    plotter->get2DHistogram("p[e-] v p[e-] - matched, time")->GetYaxis()->SetTitle("p[e-] [GeV]"); 

    plotter->build2DHistogram("p[e-] v p[e-] - moller", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e-] v p[e-] - moller")->GetXaxis()->SetTitle("p[e-] [GeV]"); 
    plotter->get2DHistogram("p[e-] v p[e-] - moller")->GetYaxis()->SetTitle("p[e-] [GeV]");
    
    // Track time
    plotter->build1DHistogram("track pair dt - matched, e-e-, time", 100, -10, 10)->GetXaxis()->SetTitle("#Delta Track time [ns]");

    plotter->build2DHistogram("track time - matched", 100, -10, 10, 100, -10, 10);
    plotter->get2DHistogram("track time - matched")->GetXaxis()->SetTitle("Track time [ns]");
    plotter->get2DHistogram("track time - matched")->GetYaxis()->SetTitle("Track time [ns]");


    plotter->build2DHistogram("track time - matched, e-e-, time", 100, -10, 10, 100, -10, 10);
    plotter->get2DHistogram("track time - matched, e-e-, time")->GetXaxis()->SetTitle("Track time [ns]");
    plotter->get2DHistogram("track time - matched, e-e-, time")->GetYaxis()->SetTitle("Track time [ns]");

    // Track theta
    plotter->build1DHistogram("delta theta - matched", 40, -0.1, 0.1); 
    
    plotter->build1DHistogram("delta theta - matched, e+e-", 40, -0.1, 0.1); 

    plotter->build2DHistogram("track theta - matched", 40, 0, 0.1, 40, 0, 0.1);
    plotter->get2DHistogram("track theta - matched")->GetXaxis()->SetTitle("#theta");
    plotter->get2DHistogram("track theta - matched")->GetYaxis()->SetTitle("#theta");


    plotter->build2DHistogram("track theta - matched, e-e-, time", 40, 0, 0.1, 40, 0, 0.1);
    plotter->get2DHistogram("track theta - matched, e-e-, time")->GetXaxis()->SetTitle("#theta");
    plotter->get2DHistogram("track theta - matched, e-e-, time")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("track theta - moller", 40, 0, 0.1, 40, 0, 0.1);
    plotter->get2DHistogram("track theta - moller")->GetXaxis()->SetTitle("#theta");
    plotter->get2DHistogram("track theta - moller")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("p v theta - matched", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - matched")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - matched")->GetYaxis()->SetTitle("#theta");


    plotter->build2DHistogram("p v theta - matched, e-e-, time", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - matched, e-e-, time")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - matched, e-e-, time")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("p v theta - moller", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - moller")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - moller")->GetYaxis()->SetTitle("#theta");

    // Invariant mass
    plotter->build1DHistogram("invariant mass - mollers", 50, 0, 0.1)->GetXaxis()->SetTitle("Invariant mass (GeV)");
    plotter->build1DHistogram("invariant mass - mollers - GBL", 50, 0, 0.1)->GetXaxis()->SetTitle("Invariant mass (GeV)");
    
    //-----------------------//
    //   Cluster positions   //
    //-----------------------//



    plotter->build2DHistogram("cluster x vs cluster x - " + cut, 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster x vs cluster x - " + cut)->GetXaxis()->SetTitle("Cluster position - x (mm)");
    plotter->get2DHistogram("cluster x vs cluster x - " + cut)->GetYaxis()->SetTitle("Cluster position - x (mm)");

    plotter->build2DHistogram("cluster x vs cluster x - moller", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster x vs cluster x - moller")->GetYaxis()->SetTitle("Cluster position - x (mm)");
    plotter->get2DHistogram("cluster x vs cluster x - moller")->GetYaxis()->SetTitle("Cluster position - x (mm)");

    plotter->build2DHistogram("cluster y vs cluster y - " + cut, 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster y vs cluster y - " + cut)->GetYaxis()->SetTitle("Cluster position - y (mm)");
    plotter->get2DHistogram("cluster y vs cluster y - " + cut)->GetYaxis()->SetTitle("Cluster position - y (mm)");

    plotter->build2DHistogram("cluster y vs cluster y - moller", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster y vs cluster y - moller")->GetYaxis()->SetTitle("Cluster position - y (mm)");
    plotter->get2DHistogram("cluster y vs cluster y - moller")->GetYaxis()->SetTitle("Cluster position - y (mm)");

    plotter->build1DHistogram("cluster pair delta x - " + cut, 50, -200, 200)->GetXaxis()->SetTitle("Ecal cluster pair dx (mm)");
    plotter->build1DHistogram("cluster pair delta x - moller", 50, -200, 200)->GetXaxis()->SetTitle("Ecal cluster pair dx (mm)");

    plotter->build1DHistogram("cluster pair x sum - " + cut, 50, -250, -100)->GetXaxis()->SetTitle("Ecal cluster pair x sum (mm)");
    plotter->build1DHistogram("cluster pair x sum - moller", 50, -250, -100)->GetXaxis()->SetTitle("Ecal cluster pair x sum (mm)");

    plotter->build2DHistogram("cluster position - moller", 200, -200, 200, 100, -100, 100);

    //------------------------------------//
    //--- Extrapolated track positions ---//
    //------------------------------------//

    // All tracks and clusters
    
    // Top
    plotter->build2DHistogram("cluster x v extrapolated track x - top", 200, -200, 200, 200, -200, 200);
    plotter->build2DHistogram("p v extrapolated track x - top", 50, 0, 1.5, 200, -200, 200);
    plotter->build2DHistogram("p v cluster x - top", 50, 0, 1.5, 200, -200, 200);
    plotter->build2DHistogram("cluster pair energy v cluster x - top", 50, 0, 1.5, 200, -200, 200);
    plotter->build2DHistogram("cluster x - track x v e/p - top", 200, -200, 200, 40, 0, 2);
   

    // Bottom
    plotter->build2DHistogram("cluster x v extrapolated track x - bottom", 200, -200, 200, 200, -200, 200);
    plotter->build2DHistogram("p v extrapolated track x - bottom", 50, 0, 1.5, 200, -200, 200);
    plotter->build2DHistogram("p v cluster x - bottom", 50, 0, 1.5, 200, -200, 200);
    plotter->build2DHistogram("cluster pair energy v cluster x - bottom", 50, 0, 1.5, 200, -200, 200);
    plotter->build2DHistogram("cluster x - track x v e/p - bottom", 200, -200, 200, 40, 0, 2);

    plotter->build2DHistogram("cluster y v extrapolated track y - bottom", 100, -100, 100, 100, -100, 100);

    // Matched tracks
    
    plotter->build2DHistogram("cluster x v e/p", 200, -200, 200, 50, 0, 1.5);
    plotter->build2DHistogram("cluster y v e/p", 50, 0, 1.5, 100, -100, 100);

    // Top
    plotter->build2DHistogram("cluster x v extrapolated track x - top - matched", 200, -200, 200, 200, -200, 200);
    plotter->build2DHistogram("cluster y v extrapolated track y - top - matched", 100, -100, 100, 100, -100, 100);
   
    // Bottom
    plotter->build2DHistogram("cluster x v extrapolated track x - bottom - matched", 200, -200, 200, 200, -200, 200);
    plotter->build2DHistogram("cluster y v extrapolated track y - bottom - matched", 100, -100, 100, 100, -100, 100);
} 

std::string MollerAnalysis::toString() { 

    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}


bool MollerAnalysis::isMatch(EcalCluster* cluster, SvtTrack* track) { 

    // Get the cluster position
    std::vector<double> cluster_pos = cluster->getPosition();

    // Extrapolate the track to the Ecal cluster position 
    std::vector<double> track_pos_at_ecal 
        = TrackExtrapolator::extrapolateTrack(track, cluster_pos[2]);

    double p = AnalysisUtils::getMagnitude(track->getMomentum());

    double delta_x = cluster_pos[0] - track_pos_at_ecal[0];
    double delta_y = cluster_pos[1] - track_pos_at_ecal[1];

    if (track->isTopTrack()) { 

        plotter->get2DHistogram("p v extrapolated track x - top")->Fill(p, track_pos_at_ecal[0]);
        plotter->get2DHistogram("p v cluster x - top")->Fill(p, cluster_pos[0]);
        plotter->get2DHistogram("cluster pair energy v cluster x - top")->Fill(cluster->getEnergy(), cluster_pos[0]);
        plotter->get2DHistogram("cluster x - track x v e/p - top")->Fill(delta_x, cluster->getEnergy()/p); 
    
    } else {
        
        plotter->get2DHistogram("p v extrapolated track x - bottom")->Fill(p, track_pos_at_ecal[0]);
        plotter->get2DHistogram("p v cluster x - bottom")->Fill(p, cluster_pos[0]);
        plotter->get2DHistogram("cluster pair energy v cluster x - bottom")->Fill(cluster->getEnergy(), cluster_pos[0]);
        plotter->get2DHistogram("cluster x - track x v e/p - bottom")->Fill(delta_x,
               cluster->getEnergy()/p); 
    }
    
    return true;
}

