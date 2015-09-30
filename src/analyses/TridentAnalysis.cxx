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
    
    // Only look at pairs1 triggers
    if (!event->isPair1Trigger()) return;

    // Only look at events with the SVT bias ON
    if (!event->isSvtBiasOn()) return; 
    
    // Only look at events where the SVT is closed
    if (!event->isSvtClosed()) return;

    // Get a "good" pair from the event.  If a good pair isn't found, skip
    // the event.
    std::vector<EcalCluster*> pair = ecal_utils->getClusterPair(event);
    if (pair.size() != 2 || pair[0] == nullptr || pair[1] == nullptr) return;
    
    // Find all track-cluster matches in the event
    matcher->findAllMatches(event);

    // Check if matches was found for the two clusters
    if (matcher->getMatchingTrack(pair[0]) == nullptr 
            || matcher->getMatchingTrack(pair[1]) == nullptr) return;

    plotter->get2DHistogram("cluster pair energy - matched")->Fill(pair[0]->getEnergy(), 
            pair[1]->getEnergy()); 
    plotter->get1DHistogram("cluster pair energy sum - matched")->Fill(
            pair[0]->getEnergy() + pair[1]->getEnergy()); 
    plotter->get1DHistogram("cluster pair energy diff - matched")->Fill(
            pair[0]->getEnergy() - pair[1]->getEnergy());

    SvtTrack* electron = nullptr;
    SvtTrack* positron = nullptr;

    if (matcher->getMatchingTrack(pair[0])->getCharge()
            *matcher->getMatchingTrack(pair[1])->getCharge() >= 0) return; 

    plotter->get2DHistogram("cluster pair energy - matched, e+e-")->Fill(pair[0]->getEnergy(), 
            pair[1]->getEnergy()); 
    plotter->get1DHistogram("cluster pair energy sum - matched, e+e-")->Fill(
            pair[0]->getEnergy() + pair[1]->getEnergy()); 
    plotter->get1DHistogram("cluster pair energy diff - matched, e+e-")->Fill(
            pair[0]->getEnergy() - pair[1]->getEnergy());


    if (matcher->getMatchingTrack(pair[0])->getCharge() == -1) { 
        electron = matcher->getMatchingTrack(pair[0]);  
        positron = matcher->getMatchingTrack(pair[1]);  
    } else { 
        electron = matcher->getMatchingTrack(pair[1]);  
        positron = matcher->getMatchingTrack(pair[0]);  
    }
    
    double electron_p = AnalysisUtils::getMagnitude(electron->getMomentum());
    double positron_p = AnalysisUtils::getMagnitude(positron->getMomentum());

    plotter->get2DHistogram("p[e+] v p[e-] - matched, e+e-")->Fill(electron_p, positron_p); 
     
    plotter->get1DHistogram("invariant mass")->Fill(AnalysisUtils::getInvariantMass(electron, positron)); 
}

void TridentAnalysis::finalize() { 
    ecal_utils->saveHistograms();
    plotter->saveToRootFile("trident_analysis.root");
}

void TridentAnalysis::bookHistograms() { 

    plotter->build2DHistogram("cluster pair energy - matched", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - matched")->GetXaxis()->SetTitle("Cluster Energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - matched")->GetYaxis()->SetTitle("Cluster Energy (GeV)");

    plotter->build2DHistogram("cluster pair energy - matched, e+e-", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - matched, e+e-")->GetXaxis()->SetTitle("Cluster Energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - matched, e+e-")->GetYaxis()->SetTitle("Cluster Energy (GeV)");

    plotter->build1DHistogram("cluster pair energy sum - matched", 50, 0, 1.5)->GetXaxis()->SetTitle(
            "Cluster Pair Energy Sum (GeV)");
    plotter->build1DHistogram("cluster pair energy sum - matched, e+e-", 50, 0, 1.5)->GetXaxis()->SetTitle(
            "Cluster Pair Energy Sum (GeV)");

    plotter->build1DHistogram("cluster pair energy diff - matched", 50, -1, 1)->GetXaxis()->SetTitle(
            "Cluster Pair Energy Difference (GeV)");
    plotter->build1DHistogram("cluster pair energy diff - matched, e+e-", 50, -1, 1)->GetXaxis()->SetTitle(
            "Cluster Pair Energy Difference (GeV)");

    // Track momentum
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
