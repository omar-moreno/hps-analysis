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
      class_name("TridentAnalysis") { 
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
    
    // Get a "good" pair from the event.  If a good pair isn't found, skip
    // the event.
    std::vector<EcalCluster*> pair = ecal_utils->getClusterPair(event);
    if (pair.size() != 2 || pair[0] == nullptr || pair[1] == nullptr) return;
    
    // Find all track-cluster matches in the event
    matcher->findAllMatches(event);

    // Check if matches for the two clusters were found
    if (matcher->getMatchingTrack(pair[0]) == nullptr 
            || matcher->getMatchingTrack(pair[1]) == nullptr) return;

    // Define some commonly used variables
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

    plotter->get2DHistogram("cluster pair energy - matched")->Fill(cluster_energy[0], cluster_energy[1]); 
    plotter->get1DHistogram("cluster pair energy sum - matched")->Fill(cluster_energy_sum); 
    plotter->get1DHistogram("cluster pair energy diff - matched")->Fill(cluster_energy_diff);
    plotter->get1DHistogram("track time - cluster time - matched")->Fill(
            matcher->getMatchingTrack(pair[0])->getTrackTime() - cluster_time[0]);  
    plotter->get1DHistogram("track time - cluster time - matched")->Fill(
            matcher->getMatchingTrack(pair[1])->getTrackTime() - cluster_time[1]);  
    plotter->get1DHistogram("track pair dt - matched")->Fill(
            matcher->getMatchingTrack(pair[0])->getTrackTime() - matcher->getMatchingTrack(pair[1])->getTrackTime()); 

    SvtTrack* electron = nullptr;
    SvtTrack* positron = nullptr;

    if (matcher->getMatchingTrack(pair[0])->getCharge()*matcher->getMatchingTrack(pair[1])->getCharge() >= 0) return; 

    if (matcher->getMatchingTrack(pair[0])->getCharge() == -1) { 
        electron = matcher->getMatchingTrack(pair[0]);  
        positron = matcher->getMatchingTrack(pair[1]);  
    } else { 
        electron = matcher->getMatchingTrack(pair[1]);  
        positron = matcher->getMatchingTrack(pair[0]);  
    }
    
    double electron_p = AnalysisUtils::getMagnitude(electron->getMomentum());
    double positron_p = AnalysisUtils::getMagnitude(positron->getMomentum());

    plotter->get2DHistogram("cluster pair energy - matched, e+e-")->Fill(cluster_energy[0], cluster_energy[1]); 
    plotter->get1DHistogram("cluster pair energy sum - matched, e+e-")->Fill(cluster_energy_sum); 
    plotter->get1DHistogram("cluster pair energy diff - matched, e+e-")->Fill(cluster_energy_diff);
    plotter->get1DHistogram("track time - cluster time - matched, e+e-")->Fill(
            matcher->getMatchingTrack(pair[0])->getTrackTime() - cluster_time[0]);  
    plotter->get1DHistogram("track time - cluster time - matched, e+e-")->Fill(
            matcher->getMatchingTrack(pair[1])->getTrackTime() - cluster_time[1]);  
    plotter->get1DHistogram("track pair dt - matched, e+e-")->Fill(
            matcher->getMatchingTrack(pair[0])->getTrackTime() - matcher->getMatchingTrack(pair[1])->getTrackTime()); 

    plotter->get2DHistogram("p[e+] v p[e-] - matched, e+e-")->Fill(electron_p, positron_p); 
     
    plotter->get1DHistogram("invariant mass")->Fill(AnalysisUtils::getInvariantMass(electron, positron)); 
}

void TridentAnalysis::finalize() { 
    ecal_utils->saveHistograms();
    plotter->saveToRootFile("trident_analysis.root");
}

void TridentAnalysis::bookHistograms() { 

    // Matched Ecal clusters
    plotter->build1DHistogram("cluster pair energy diff - matched", 50, -1, 1)->GetXaxis()->SetTitle(
            "Cluster Pair Energy Difference (GeV)");
    
    plotter->build1DHistogram("track time - cluster time - matched", 60, 30, 60)->GetXaxis()->SetTitle(
            "Track time - Ecal cluster time"); 

    plotter->build2DHistogram("cluster pair energy - matched", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - matched")->GetXaxis()->SetTitle("Cluster Energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - matched")->GetYaxis()->SetTitle("Cluster Energy (GeV)");

    plotter->build1DHistogram("cluster pair energy sum - matched", 50, 0, 1.5)->GetXaxis()->SetTitle(
            "Cluster Pair Energy Sum (GeV)");
    
    plotter->build1DHistogram("cluster pair energy diff - matched, e+e-", 50, -1, 1)->GetXaxis()->SetTitle(
            "Cluster Pair Energy Difference (GeV)");

    plotter->build1DHistogram("track time - cluster time - matched, e+e-", 60, 30, 60)->GetXaxis()->SetTitle(
            "Track time - Ecal cluster time");

    plotter->build1DHistogram("cluster pair energy sum - matched, e+e-", 50, 0, 1.5)->GetXaxis()->SetTitle(
            "Cluster Pair Energy Sum (GeV)");

    plotter->build2DHistogram("cluster pair energy - matched, e+e-", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - matched, e+e-")->GetXaxis()->SetTitle("Cluster Energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - matched, e+e-")->GetYaxis()->SetTitle("Cluster Energy (GeV)");

    // Matched tracks
    plotter->build1DHistogram("track pair dt - matched", 50, -25, 25)->GetXaxis()->SetTitle("Track time #Delta t");  
    plotter->build1DHistogram("track pair dt - matched, e+e-", 50, -25, 25)->GetXaxis()->SetTitle(
            "Track time #Delta t");  
    plotter->build1DHistogram("p top - matched, e+e-", 50, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");
    plotter->build1DHistogram("p bottom - matched, e+e-", 50, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");

    plotter->build2DHistogram("p[e+] v p[e-] - matched, e+e-", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e+] v p[e-] - matched, e+e-")->GetXaxis()->SetTitle("p[e-] (GeV)"); 
    plotter->get2DHistogram("p[e+] v p[e-] - matched, e+e-")->GetYaxis()->SetTitle("p[e+] (GeV)"); 

    // Invariant mass
    plotter->build1DHistogram("invariant mass", 50, 0, 0.1)->GetXaxis()->SetTitle("Mass (GeV)");

}

std::string TridentAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
