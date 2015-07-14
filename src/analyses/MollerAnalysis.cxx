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
      class_name("MollerAnalysis") {

}

MollerAnalysis::~MollerAnalysis() {
    delete plotter;
}

void MollerAnalysis::initialize() { 
    this->bookHistograms(); 
}

void MollerAnalysis::processEvent(HpsEvent* event) { 

    // Only look at pair1 triggers
    if (!event->isPair1Trigger()) return;

    // Only look at events that have two tracks 
    if (event->getNumberOfTracks() != 2) return;

    matcher->findAllMatches(event); 

    SvtTrack* first_track = event->getTrack(0);
    SvtTrack* second_track = event->getTrack(1); 

    // Require that both tracks are negatively charged
    if ((first_track->getCharge() + second_track->getCharge()) != -2) return;
       
    double p0 = AnalysisUtils::getMagnitude(first_track->getMomentum());
    double p1 = AnalysisUtils::getMagnitude(second_track->getMomentum()); 
    
    if (first_track->isTopTrack()) { 
        plotter->get1DHistogram("p top")->Fill(p0);
    } else { 
        plotter->get1DHistogram("p bottom")->Fill(p0);
    }

    if (second_track->isTopTrack()) { 
        plotter->get1DHistogram("p top")->Fill(p0);
    } else { 
        plotter->get1DHistogram("p bottom")->Fill(p0);
    }

    plotter->get2DHistogram("p[e-] v p[e-] - all")->Fill(p0, p1);
    
    plotter->get1DHistogram("delta Track Time - all")->Fill(first_track->getTrackTime() - second_track->getTrackTime());

    plotter->get1DHistogram("delta cos(theta)")->Fill(TrackExtrapolator::getCosTheta(first_track) - 
            TrackExtrapolator::getCosTheta(second_track));

    // Require that the tracks are in opposite volumes
    if ((first_track->isTopTrack() && second_track->isTopTrack()) 
            || (first_track->isBottomTrack() && second_track->isBottomTrack())) return;

    // Require that both tracks are matched to a cluster
    if (matcher->getMatchingCluster(first_track) == NULL || 
            matcher->getMatchingCluster(second_track) == NULL) return;
    plotter->get2DHistogram("p[e-] v p[e-] - matched")->Fill(p0, p1);
    plotter->get1DHistogram("delta Track Time - matched")->Fill(first_track->getTrackTime() - second_track->getTrackTime());

    EcalCluster* first_cluster = matcher->getMatchingCluster(first_track);
    EcalCluster* second_cluster = matcher->getMatchingCluster(second_track);


    plotter->get2DHistogram("first v second cluster time")->Fill(
            first_cluster->getClusterTime(), 
            second_cluster->getClusterTime());

    double delta_cluster_time = first_cluster->getClusterTime() 
        - second_cluster->getClusterTime();

    plotter->get1DHistogram("delta clusters - matched")->Fill(delta_cluster_time);

    // Require that the sum of the momentum be above and below 
    // some reasonable value
    if ((p0+p1) < .80*1.056) return;
    if ((p0+p1) > 1.2*1.056) return;

    //if (abs(first_track->getTrackTime() - second_track->getTrackTime()) > 4) return;

    //if (TrackExtrapolator::getCosTheta(first_track) - 
    //        TrackExtrapolator::getCosTheta(second_track) < 0.04) return;

    plotter->get2DHistogram("p[e-] v p[e-] - Pass Cuts")->Fill(p0, p1);
    plotter->get1DHistogram("delta Track Time - pass cuts")->Fill(first_track->getTrackTime() - second_track->getTrackTime());
    plotter->get1DHistogram("delta clusters - pass cuts")->Fill(delta_cluster_time);
    
}

void MollerAnalysis::finalize() { 
    plotter->saveToPdf("moller_analysis.pdf");
}


void MollerAnalysis::bookHistograms() {

    plotter->build1DHistogram("p top", 50, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");
    plotter->build1DHistogram("p bottom", 50, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");

    plotter->build1DHistogram("delta Track Time - all", 100, -10, 10)->GetXaxis()->SetTitle("#Delta Track time [ns]");
    plotter->build1DHistogram("delta Track Time - matched", 100, -10, 10)->GetXaxis()->SetTitle("#Delta Track time [ns]");
    plotter->build1DHistogram("delta Track Time - pass cuts", 100, -10, 10)->GetXaxis()->SetTitle("#Delta Track time [ns]");

    plotter->build1DHistogram("delta clusters - matched", 100, -10, 10)->GetXaxis()->SetTitle("#Delta Cluster time [ns]");

    plotter->build1DHistogram("delta clusters - pass cuts", 100, -10, 10)->GetXaxis()->SetTitle("#Delta Cluster time [ns]");

    plotter->build1DHistogram("delta cos(theta)", 40, -0.1, 0.1); 

    plotter->build2DHistogram("p[e-] v p[e-] - all", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e-] v p[e-] - all")->GetXaxis()->SetTitle("p[e-] [GeV]"); 
    plotter->get2DHistogram("p[e-] v p[e-] - all")->GetYaxis()->SetTitle("p[e-] [GeV]"); 

    plotter->build2DHistogram("p[e-] v p[e-] - matched", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e-] v p[e-] - matched")->GetXaxis()->SetTitle("p[e-] [GeV]"); 
    plotter->get2DHistogram("p[e-] v p[e-] - matched")->GetYaxis()->SetTitle("p[e-] [GeV]"); 

    plotter->build2DHistogram("first v second cluster time", 160, 20, 100, 160, 20, 100);
    plotter->get2DHistogram("first v second cluster time")->GetXaxis()->SetTitle("First cluster time (ns)");
    plotter->get2DHistogram("first v second cluster time")->GetYaxis()->SetTitle("Second cluster time (ns)");

    plotter->build2DHistogram("first v second cluster energy", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("first v second cluster energy")->GetYaxis()->SetTitle("First cluster energy (GeV)");
    plotter->get2DHistogram("first v second cluster energy")->GetYaxis()->SetTitle("Second cluster energy (GeV)");

    plotter->build2DHistogram("p[e-] v p[e-] - Pass Cuts", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e-] v p[e-] - Pass Cuts")->GetXaxis()->SetTitle("p[e-] [GeV]"); 
    plotter->get2DHistogram("p[e-] v p[e-] - Pass Cuts")->GetYaxis()->SetTitle("p[e-] [GeV]"); 
} 

std::string MollerAnalysis::toString() { 
    
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
