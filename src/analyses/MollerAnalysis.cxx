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
      class_name("MollerAnalysis") {

}

MollerAnalysis::~MollerAnalysis() {
    delete plotter;
}

void MollerAnalysis::initialize() { 
    this->bookHistograms(); 
}

void MollerAnalysis::processEvent(HpsEvent* event) { 

    // Only look at events that have two tracks 
    if (event->getNumberOfTracks() != 2) return;

    SvtTrack* first_track = event->getTrack(0);
    SvtTrack* second_track = event->getTrack(1); 

    if (first_track->getCharge()*second_track->getCharge() != 1) return;
       
    double p0 = AnalysisUtils::getMagnitude(first_track->getMomentum());
    double p1 = AnalysisUtils::getMagnitude(second_track->getMomentum()); 
        
    plotter->get2DHistogram("p[e-] v p[e-]")->Fill(p0, p1);
    
    plotter->get1DHistogram("delta Track Time")->Fill(first_track->getTrackTime() - second_track->getTrackTime());

    plotter->get1DHistogram("delta cos(theta)")->Fill(TrackExtrapolator::getCosTheta(first_track) - 
            TrackExtrapolator::getCosTheta(second_track));

    if ((first_track->isTopTrack() && second_track->isTopTrack()) 
            || (first_track->isBottomTrack() && second_track->isBottomTrack())) return;

    if (p0 > .8*1.056 || p1 > .8*1.056) return;

    if (abs(first_track->getTrackTime() - second_track->getTrackTime()) > 4) return;

    if (TrackExtrapolator::getCosTheta(first_track) - 
            TrackExtrapolator::getCosTheta(second_track) < 0.04) return;

    plotter->get2DHistogram("p[e-] v p[e-] - Pass Cuts")->Fill(p0, p1);
    
}

void MollerAnalysis::finalize() { 
    plotter->saveToPdf("moller_analysis.pdf");
}


void MollerAnalysis::bookHistograms() {

    plotter->build1DHistogram("delta Track Time", 100, -10, 10)->GetXaxis()->SetTitle("delta Track time [ns]");
    plotter->build1DHistogram("delta cos(theta)", 40, -0.1, 0.1); 

    plotter->build2DHistogram("p[e-] v p[e-]", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e-] v p[e-]")->GetXaxis()->SetTitle("p[e-] [GeV]"); 
    plotter->get2DHistogram("p[e-] v p[e-]")->GetYaxis()->SetTitle("p[e-] [GeV]"); 

    plotter->build2DHistogram("p[e-] v p[e-] - Pass Cuts", 50, 0, 1., 50, 0, 1.);
    plotter->get2DHistogram("p[e-] v p[e-] - Pass Cuts")->GetXaxis()->SetTitle("p[e-] [GeV]"); 
    plotter->get2DHistogram("p[e-] v p[e-] - Pass Cuts")->GetYaxis()->SetTitle("p[e-] [GeV]"); 
} 

std::string MollerAnalysis::toString() { 
    
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
