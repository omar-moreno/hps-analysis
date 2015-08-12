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

        // Calculate the momentum magnitude and transverse momentum
        std::vector<double> gbl_p = gbl_track->getMomentum();
        double gbl_p_mag = sqrt(gbl_p[0]*gbl_p[0] + gbl_p[1]*gbl_p[1] + gbl_p[2]*gbl_p[2]);
        double gbl_pt = sqrt(gbl_p[0]*gbl_p[0] + gbl_p[1]*gbl_p[1]);

        plotter->get1DHistogram("p - gbl")->Fill(gbl_p_mag);
        plotter->get1DHistogram("pt - gbl")->Fill(gbl_pt);
        plotter->get1DHistogram("px - gbl")->Fill(gbl_p[0]);
        plotter->get1DHistogram("py - gbl")->Fill(gbl_p[1]);
        plotter->get1DHistogram("pz - gbl")->Fill(gbl_p[2]);

        plotter->get1DHistogram("chi2")->Fill(seed_track->getChi2());
        plotter->get1DHistogram("doca")->Fill(seed_track->getD0());
        plotter->get1DHistogram("z0")->Fill(seed_track->getZ0());
        plotter->get1DHistogram("sin(phi0)")->Fill(sin(seed_track->getPhi0()));
        plotter->get1DHistogram("curvature")->Fill(seed_track->getOmega());
        plotter->get1DHistogram("cos(theta)")->Fill(TrackExtrapolator::getCosTheta(seed_track));

        // Calculate the momentum magnitude and transverse momentum
        std::vector<double> p = seed_track->getMomentum();
        double p_mag = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        double pt = sqrt(p[0]*p[0] + p[1]*p[1]);

        plotter->get1DHistogram("p")->Fill(p_mag);
        plotter->get1DHistogram("pt")->Fill(pt);
        plotter->get1DHistogram("px")->Fill(p[0]);
        plotter->get1DHistogram("py")->Fill(p[1]);
        plotter->get1DHistogram("pz")->Fill(p[2]);


        if (seed_track->isTopTrack()) { 
            
            plotter->get1DHistogram("chi2 - top")->Fill(seed_track->getChi2()); 
            plotter->get1DHistogram("doca - top")->Fill(seed_track->getD0());
            plotter->get1DHistogram("z0 - top")->Fill(seed_track->getZ0());
            plotter->get1DHistogram("sin(phi0) - top")->Fill(sin(seed_track->getPhi0()));
            plotter->get1DHistogram("curvature - top")->Fill(seed_track->getOmega());
            plotter->get1DHistogram("cos(theta) - top")->Fill(TrackExtrapolator::getCosTheta(seed_track));
        
            plotter->get1DHistogram("p - top")->Fill(p_mag);
            plotter->get1DHistogram("pt - top")->Fill(pt);
            plotter->get1DHistogram("px - top")->Fill(p[0]);
            plotter->get1DHistogram("py - top")->Fill(p[1]);
            plotter->get1DHistogram("pz - top")->Fill(p[2]);

            plotter->get1DHistogram("chi2 - top - gbl")->Fill(gbl_track->getChi2()); 
            plotter->get1DHistogram("doca - top - gbl")->Fill(gbl_track->getD0());
            plotter->get1DHistogram("z0 - top - gbl")->Fill(gbl_track->getZ0());
            plotter->get1DHistogram("sin(phi0) - top - gbl")->Fill(sin(gbl_track->getPhi0()));
            plotter->get1DHistogram("curvature - top - gbl")->Fill(gbl_track->getKappa());
            plotter->get1DHistogram("cos(theta) - top - gbl")->Fill(cos(gbl_track->getTheta()));
        
            plotter->get1DHistogram("p - top - gbl")->Fill(gbl_p_mag);
            plotter->get1DHistogram("pt - top - gbl")->Fill(gbl_pt);
            plotter->get1DHistogram("px - top - gbl")->Fill(gbl_p[0]);
            plotter->get1DHistogram("py - top - gbl")->Fill(gbl_p[1]);
            plotter->get1DHistogram("pz - top - gbl")->Fill(gbl_p[2]);
        
        } else { 
            
            plotter->get1DHistogram("chi2 - bottom")->Fill(seed_track->getChi2()); 
            plotter->get1DHistogram("doca - bottom")->Fill(seed_track->getD0());
            plotter->get1DHistogram("z0 - bottom")->Fill(seed_track->getZ0());
            plotter->get1DHistogram("sin(phi0) - bottom")->Fill(sin(seed_track->getPhi0()));
            plotter->get1DHistogram("curvature - bottom")->Fill(seed_track->getOmega());
            plotter->get1DHistogram("cos(theta) - bottom")->Fill(TrackExtrapolator::getCosTheta(seed_track));
        
            plotter->get1DHistogram("p - bottom")->Fill(p_mag);
            plotter->get1DHistogram("pt - bottom")->Fill(pt);
            plotter->get1DHistogram("px - bottom")->Fill(p[0]);
            plotter->get1DHistogram("py - bottom")->Fill(p[1]);
            plotter->get1DHistogram("pz - bottom")->Fill(p[2]);

            plotter->get1DHistogram("chi2 - bottom - gbl")->Fill(gbl_track->getChi2()); 
            plotter->get1DHistogram("doca - bottom - gbl")->Fill(gbl_track->getD0());
            plotter->get1DHistogram("z0 - bottom - gbl")->Fill(gbl_track->getZ0());
            plotter->get1DHistogram("sin(phi0) - bottom - gbl")->Fill(sin(gbl_track->getPhi0()));
            plotter->get1DHistogram("curvature - bottom - gbl")->Fill(gbl_track->getKappa());
            plotter->get1DHistogram("cos(theta) - bottom - gbl")->Fill(cos(gbl_track->getTheta()));
        
            plotter->get1DHistogram("p - bottom - gbl")->Fill(gbl_p_mag);
            plotter->get1DHistogram("pt - bottom - gbl")->Fill(gbl_pt);
            plotter->get1DHistogram("px - bottom - gbl")->Fill(gbl_p[0]);
            plotter->get1DHistogram("py - bottom - gbl")->Fill(gbl_p[1]);
            plotter->get1DHistogram("pz - bottom - gbl")->Fill(gbl_p[2]);
        
        }


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

    //-------------------//
    //--- Seed Tracks ---//
    //-------------------//

    plotter->build1DHistogram("chi2", 25, 0, 25)->GetXaxis()->SetTitle("Track #chi^{2}");
    plotter->build1DHistogram("doca", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("z0", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("sin(phi0)", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("curvature", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("cos(theta)", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    plotter->build1DHistogram("p", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    plotter->build1DHistogram("pt", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{t} [GeV]");
    plotter->build1DHistogram("px", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{x} [GeV]"); 
    plotter->build1DHistogram("py", 50, -0.15, 0.15)->GetXaxis()->SetTitle("p_{y} [GeV]"); 
    plotter->build1DHistogram("pz", 50, 0, 2.0)->GetXaxis()->SetTitle("p_{z} [GeV]");

    // Top tracks
    plotter->build1DHistogram("chi2 - top", 25, 0, 25)->GetXaxis()->SetTitle("Track #chi^{2}");
    plotter->build1DHistogram("doca - top", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("z0 - top", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("sin(phi0) - top", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("curvature - top", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("cos(theta) - top", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    plotter->build1DHistogram("p - top", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    plotter->build1DHistogram("pt - top", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{t} [GeV]");
    plotter->build1DHistogram("px - top", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{x} [GeV]"); 
    plotter->build1DHistogram("py - top", 50, -0.15, 0.15)->GetXaxis()->SetTitle("p_{y} [GeV]"); 
    plotter->build1DHistogram("pz - top", 50, 0, 2.0)->GetXaxis()->SetTitle("p_{z} [GeV]");

    // Bottom tracks
    plotter->build1DHistogram("chi2 - bottom", 25, 0, 25)->GetXaxis()->SetTitle("Track #chi^{2}");
    plotter->build1DHistogram("doca - bottom", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("z0 - bottom", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("sin(phi0) - bottom", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("curvature - bottom", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("cos(theta) - bottom", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    plotter->build1DHistogram("p - bottom", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    plotter->build1DHistogram("pt - bottom", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{t} [GeV]");
    plotter->build1DHistogram("px - bottom", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{x} [GeV]"); 
    plotter->build1DHistogram("py - bottom", 50, -0.15, 0.15)->GetXaxis()->SetTitle("p_{y} [GeV]"); 
    plotter->build1DHistogram("pz - bottom", 50, 0, 2.0)->GetXaxis()->SetTitle("p_{z} [GeV]");

    //------------------//
    //--- GBL Tracks ---//
    //------------------//

    plotter->build1DHistogram("p - gbl", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    plotter->build1DHistogram("pt - gbl", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{t} [GeV]");
    plotter->build1DHistogram("px - gbl", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{x} [GeV]"); 
    plotter->build1DHistogram("py - gbl", 50, -0.15, 0.15)->GetXaxis()->SetTitle("p_{y} [GeV]"); 
    plotter->build1DHistogram("pz - gbl", 50, 0, 2.0)->GetXaxis()->SetTitle("p_{z} [GeV]");

    plotter->build1DHistogram("chi2 - gbl", 25, 0, 25)->GetXaxis()->SetTitle("Track #chi^{2}");
    plotter->build1DHistogram("doca - gbl", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("z0 - gbl", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("sin(phi0) - gbl", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("curvature - gbl", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("cos(theta) - gbl", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    // Top tracks
    plotter->build1DHistogram("chi2 - top - gbl", 25, 0, 25)->GetXaxis()->SetTitle("Track #chi^{2}");
    plotter->build1DHistogram("doca - top - gbl", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("z0 - top - gbl", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("sin(phi0) - top - gbl", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("curvature - top - gbl", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("cos(theta) - top - gbl", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    plotter->build1DHistogram("p - top - gbl", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    plotter->build1DHistogram("pt - top - gbl", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{t} [GeV]");
    plotter->build1DHistogram("px - top - gbl", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{x} [GeV]"); 
    plotter->build1DHistogram("py - top - gbl", 50, -0.15, 0.15)->GetXaxis()->SetTitle("p_{y} [GeV]"); 
    plotter->build1DHistogram("pz - top - gbl", 50, 0, 2.0)->GetXaxis()->SetTitle("p_{z} [GeV]");

    // Bottom tracks
    plotter->build1DHistogram("chi2 - bottom - gbl", 25, 0, 25)->GetXaxis()->SetTitle("Track #chi^{2}");
    plotter->build1DHistogram("doca - bottom - gbl", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("z0 - bottom - gbl", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("sin(phi0) - bottom - gbl", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("curvature - bottom - gbl", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("cos(theta) - bottom - gbl", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    plotter->build1DHistogram("p - bottom - gbl", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    plotter->build1DHistogram("pt - bottom - gbl", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{t} [GeV]");
    plotter->build1DHistogram("px - bottom - gbl", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{x} [GeV]"); 
    plotter->build1DHistogram("py - bottom - gbl", 50, -0.15, 0.15)->GetXaxis()->SetTitle("p_{y} [GeV]"); 
    plotter->build1DHistogram("pz - bottom - gbl", 50, 0, 2.0)->GetXaxis()->SetTitle("p_{z} [GeV]");

    //-------------------//
    //--- Comparisons ---//
    //-------------------//

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
