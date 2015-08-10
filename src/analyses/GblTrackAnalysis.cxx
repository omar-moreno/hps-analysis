/**
 * @file GblTrackAnalysis.cxx
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date August 4, 2015
 */

#include <GblTrackAnalysis.h>



GblTrackAnalysis::GblTrackAnalysis()
    : plotter(new Plotter()),
      class_name("GblTrackAnalysis") {  
}

GblTrackAnalysis::~GblTrackAnalysis() { 
    delete plotter; 
}

void GblTrackAnalysis::initialize() { 
    this->bookHistograms();
}

void GblTrackAnalysis::processEvent(HpsEvent* event) {

    if (!event->isSingle1Trigger()) return;

    // Loop over all of the tracks in the event
    for (int gbl_track_n = 0; gbl_track_n < event->getNumberOfGblTracks(); ++gbl_track_n) { 

        // Get a GBL track from the event
        GblTrack* gbl_track = event->getGblTrack(gbl_track_n); 
        //std::cout << "GBL Track: D0: " << gbl_track->getD0() << std::endl;
        
        // Get the seed track associated with the GBL track
        SvtTrack* seed_track = (SvtTrack*) gbl_track->getSeedTrack().GetObject();
        //std::cout << "Seed Track: D0: " << seed_track->getD0() << std::endl;

        // Fill the plots providing general event information related to tracks
        //int track_volume = track->isTopTrack() ? 0 : 1; 
       
        plotter->get1DHistogram("chi2 - gbl")->Fill(gbl_track->getChi2()); 
        plotter->get1DHistogram("doca - gbl")->Fill(gbl_track->getD0());
        plotter->get1DHistogram("z0 - gbl")->Fill(gbl_track->getZ0());
        plotter->get1DHistogram("sin(phi0) - gbl")->Fill(sin(gbl_track->getPhi0()));
        plotter->get1DHistogram("curvature - gbl")->Fill(gbl_track->getKappa());
        plotter->get1DHistogram("cos(theta) - gbl")->Fill(cos(gbl_track->getTheta()));
    
        plotter->get1DHistogram("chi2")->Fill(seed_track->getChi2());
        plotter->get1DHistogram("doca")->Fill(seed_track->getD0());
        plotter->get1DHistogram("z0")->Fill(seed_track->getZ0());
        plotter->get1DHistogram("sin(phi0)")->Fill(sin(seed_track->getPhi0()));
        plotter->get1DHistogram("curvature")->Fill(seed_track->getOmega());
        plotter->get1DHistogram("cos(theta)")->Fill(TrackExtrapolator::getCosTheta(seed_track));

        plotter->get1DHistogram("seed doca - gbl doca")->Fill(
                seed_track->getD0() - gbl_track->getD0());
        plotter->get1DHistogram("seed z0 - gbl z0")->Fill(
                seed_track->getZ0() - gbl_track->getZ0()); 
        plotter->get1DHistogram("seed sin(phi0) - gbl sin(phi0)")->Fill(
                sin(seed_track->getPhi0()) - sin(gbl_track->getPhi0())); 
        plotter->get1DHistogram("seed curvature - gbl curvature")->Fill(
                seed_track->getOmega() - gbl_track->getKappa()); 
        plotter->get1DHistogram("seed cos(theta) - gbl cos(theta)")->Fill(
                TrackExtrapolator::getCosTheta(seed_track) - cos(gbl_track->getTheta()));

    } 
}

void GblTrackAnalysis::finalize() {
    
    plotter->saveToPdf("gbl_analysis.pdf");
    plotter->saveToRootFile("gbl_analysis.root");
}

void GblTrackAnalysis::bookHistograms() {

    plotter->build1DHistogram("chi2", 25, 0, 25)->GetXaxis()->SetTitle("Track #chi^{2}");
    plotter->build1DHistogram("doca", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("z0", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("sin(phi0)", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("curvature", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("cos(theta)", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    plotter->build1DHistogram("chi2 - gbl", 25, 0, 25)->GetXaxis()->SetTitle("Track #chi^{2}");
    plotter->build1DHistogram("doca - gbl", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("z0 - gbl", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("sin(phi0) - gbl", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("curvature - gbl", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("cos(theta) - gbl", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    plotter->build1DHistogram("seed chi2 - gbl chi2", 100, -3, 3)->GetXaxis()->SetTitle("Seed #chi^{2} - GBL #chi^{2}");
    plotter->build1DHistogram("seed doca - gbl doca", 100, -2, 2)->GetXaxis()->SetTitle("Seed d0 - GBL d0");
    plotter->build1DHistogram("seed z0 - gbl z0", 100, -1, 1)->GetXaxis()->SetTitle("Seed z0 - GBL z0");
    plotter->build1DHistogram("seed sin(phi0) - gbl sin(phi0)", 100, -0.0001, 0.0001)->GetXaxis()->SetTitle("Seed sin(#phi_{0}) - GBL sin(#phi_{0})");
    plotter->build1DHistogram("seed curvature - gbl curvature", 100, -.0001, .0001)->GetXaxis()->SetTitle("Seed #Omega - GBL #Omega");
    plotter->build1DHistogram("seed cos(theta) - gbl cos(theta)", 20, -5, 5)->GetXaxis()->SetTitle("Seed cos(#theta) - GBL cos(#theta)");
} 

std::string GblTrackAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;   
}
