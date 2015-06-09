/**
 * @file TrackAnalysis.cxx
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date May 14, 2015
 *
 */

#include <TrackAnalysis.h>

TrackAnalysis::TrackAnalysis()
    : track(NULL),
      track_plotter(new Plotter()),
      electron_plotter(new Plotter()),
      positron_plotter(new Plotter()),
      top_plotter(new Plotter()), 
      bottom_plotter(new Plotter()), 
      class_name("TrackAnalysis") {  
}

TrackAnalysis::~TrackAnalysis() { 
    delete track_plotter; 
    delete electron_plotter;
    delete positron_plotter; 
    delete top_plotter;
    delete bottom_plotter;
}

void TrackAnalysis::initialize() { 

    this->bookHistograms();
}

void TrackAnalysis::processEvent(HpsEvent* event) {

    if (!event->isPair1Trigger()) return;

    track_plotter->get1DHistogram("Number of tracks")->Fill(event->getNumberOfTracks());

    // Loop over all of the tracks in the event
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 

        // Get a track from the event
        track = event->getTrack(track_n); 

        //if (track->getD0() < 2) continue;
        //if (track->getChi2() > 1) continue;

        // Fill the plots providing general event information related to tracks
        int track_volume = track->isTopTrack() ? 0 : 1; 
        track_plotter->get1DHistogram("Track volume")->Fill(track_volume);
        track_plotter->get1DHistogram("Track charge")->Fill(track->getCharge());
        track_plotter->get1DHistogram("chi2")->Fill(track->getChi2());

        // Fill the track parameter plots for all tracks
        track_plotter->get1DHistogram("doca")->Fill(track->getD0());
        track_plotter->get1DHistogram("z0")->Fill(track->getZ0());
        track_plotter->get1DHistogram("sin(phi0)")->Fill(sin(track->getPhi0()));
        track_plotter->get1DHistogram("curvature")->Fill(track->getOmega());
        track_plotter->get1DHistogram("tan_lambda")->Fill(track->getTanLambda());
       
         
        // Calculate the momentum magnitude and transverse momentum
        std::vector<double> p = track->getMomentum();
        double p_mag = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        double pt = sqrt(p[0]*p[0] + p[1]*p[1]);
        
        // Fill the momentum plots
        track_plotter->get1DHistogram("p")->Fill(p_mag);
        track_plotter->get1DHistogram("pt")->Fill(pt);
        track_plotter->get1DHistogram("px")->Fill(p[0]);
        track_plotter->get1DHistogram("py")->Fill(p[1]);
        track_plotter->get1DHistogram("pz")->Fill(p[2]);
        track_plotter->get2DHistogram("pz v px")->Fill(p[2], p[0]);
        track_plotter->get2DHistogram("pz v py")->Fill(p[2], p[1]);
        track_plotter->get2DHistogram("p v pt")->Fill(p_mag, pt);

        // Extrapolate the track to the target
        std::vector<double> track_pos_target = TrackExtrapolator::extrapolateTrack(track, 0);
        
        track_plotter->get2DHistogram("track position at target")->Fill(track_pos_target[0], track_pos_target[1]);

        // Fill the electron and positron plots
        if (track->getCharge() < 0) { 
            electron_plotter->get1DHistogram("p")->Fill(p_mag);
            electron_plotter->get1DHistogram("pt")->Fill(pt);
            electron_plotter->get1DHistogram("px")->Fill(p[0]);
            electron_plotter->get1DHistogram("py")->Fill(p[1]);
            electron_plotter->get1DHistogram("pz")->Fill(p[2]);
            electron_plotter->get1DHistogram("chi2")->Fill(track->getChi2());

            // Fill the electron track parameter plots
            electron_plotter->get1DHistogram("doca")->Fill(track->getD0());
            electron_plotter->get1DHistogram("z0")->Fill(track->getZ0());
            electron_plotter->get1DHistogram("sin(phi0)")->Fill(sin(track->getPhi0()));
            electron_plotter->get1DHistogram("curvature")->Fill(track->getOmega());
            electron_plotter->get1DHistogram("tan_lambda")->Fill(track->getTanLambda());

            if (((HpsParticle*) track->getParticle().GetObject())->getClusters()->GetEntriesFast() != 0) { 
               double energy = ((EcalCluster*) ((HpsParticle*) track->getParticle().GetObject())->getClusters()->At(0))->getEnergy();
                electron_plotter->get1DHistogram("ep")->Fill(energy/p_mag);
            }

            electron_plotter->get2DHistogram("track position at target")->Fill(track_pos_target[0], track_pos_target[1]);

        } else { 
            positron_plotter->get1DHistogram("p")->Fill(p_mag);
            positron_plotter->get1DHistogram("pt")->Fill(pt);
            positron_plotter->get1DHistogram("px")->Fill(p[0]);
            positron_plotter->get1DHistogram("py")->Fill(p[1]);
            positron_plotter->get1DHistogram("pz")->Fill(p[2]);
            positron_plotter->get1DHistogram("chi2")->Fill(track->getChi2());
            positron_plotter->get1DHistogram("doca")->Fill(track->getD0());
            positron_plotter->get1DHistogram("z0")->Fill(track->getZ0());
            positron_plotter->get1DHistogram("sin(phi0)")->Fill(sin(track->getPhi0()));
            positron_plotter->get1DHistogram("curvature")->Fill(track->getOmega());
            positron_plotter->get1DHistogram("tan_lambda")->Fill(track->getTanLambda());

            positron_plotter->get2DHistogram("track position at target")->Fill(track_pos_target[0], track_pos_target[1]);
        }
    
        // Fill the top and bottom track plots
        if (track->isTopTrack()) { 
            
            top_plotter->get1DHistogram("p")->Fill(p_mag);
            top_plotter->get1DHistogram("pt")->Fill(pt);
            top_plotter->get1DHistogram("px")->Fill(p[0]);
            top_plotter->get1DHistogram("py")->Fill(p[1]);
            top_plotter->get1DHistogram("pz")->Fill(p[2]);
            top_plotter->get1DHistogram("chi2")->Fill(track->getChi2());
            top_plotter->get1DHistogram("doca")->Fill(track->getD0());
            top_plotter->get1DHistogram("z0")->Fill(track->getZ0());
            top_plotter->get1DHistogram("sin(phi0)")->Fill(sin(track->getPhi0()));
            top_plotter->get1DHistogram("curvature")->Fill(track->getOmega());
            top_plotter->get1DHistogram("tan_lambda")->Fill(track->getTanLambda());

            if (track->getCharge() < 0 && ((HpsParticle*) track->getParticle().GetObject())->getClusters()->GetEntriesFast() != 0) { 
               double energy = ((EcalCluster*) ((HpsParticle*) track->getParticle().GetObject())->getClusters()->At(0))->getEnergy();
                top_plotter->get1DHistogram("ep")->Fill(energy/p_mag);
            }
        
        } else { 
            
            bottom_plotter->get1DHistogram("p")->Fill(p_mag);
            bottom_plotter->get1DHistogram("pt")->Fill(pt);
            bottom_plotter->get1DHistogram("px")->Fill(p[0]);
            bottom_plotter->get1DHistogram("py")->Fill(p[1]);
            bottom_plotter->get1DHistogram("pz")->Fill(p[2]);
            bottom_plotter->get1DHistogram("chi2")->Fill(track->getChi2());
            bottom_plotter->get1DHistogram("doca")->Fill(track->getD0());
            bottom_plotter->get1DHistogram("z0")->Fill(track->getZ0());
            bottom_plotter->get1DHistogram("sin(phi0)")->Fill(sin(track->getPhi0()));
            bottom_plotter->get1DHistogram("curvature")->Fill(track->getOmega());
            bottom_plotter->get1DHistogram("tan_lambda")->Fill(track->getTanLambda());

            if (track->getCharge() < 0 && ((HpsParticle*) track->getParticle().GetObject())->getClusters()->GetEntriesFast() != 0) { 
               double energy = ((EcalCluster*) ((HpsParticle*) track->getParticle().GetObject())->getClusters()->At(0))->getEnergy();
                bottom_plotter->get1DHistogram("ep")->Fill(energy/p_mag);
            }
        }

        track_plotter->get1DHistogram("Hits per track")->Fill(track->getSvtHits()->GetEntriesFast()); 
    
        TRefArray* hits = track->getSvtHits();
        for (int hit_n = 0; hit_n < hits->GetSize(); ++hit_n) { 
            std::vector<double> position = ((SvtHit*) hits->At(hit_n))->getPosition();

            // Extrapolate the track to the target
            std::vector<double> track_pos_layer = TrackExtrapolator::extrapolateTrack(track, position[0]);

            double x_residual = track_pos_layer[0] - position[1]; 
            double y_residual = track_pos_layer[1] - position[2];

            int layer = ((SvtHit*) hits->At(hit_n))->getLayer();
            if (track->isTopTrack()) { 
                track_plotter->get2DHistogram("Top Layer " + std::to_string(layer) + " - Hit Positions")->Fill(
                        position[1], position[2]);
                track_plotter->get1DHistogram("Top Layer " + std::to_string(layer) + " - x Residuals")->Fill(
                        x_residual); 
                track_plotter->get1DHistogram("Top Layer " + std::to_string(layer) + " - y Residuals")->Fill(
                        y_residual); 
                if (track->getCharge() < 0) { 
                    electron_plotter->get2DHistogram("Top Layer " + std::to_string(layer) + " - Hit Positions")->Fill(
                        position[1], position[2]);
                    electron_plotter->get1DHistogram("Top Layer " + std::to_string(layer) + " - x Residuals")->Fill(
                        x_residual); 
                    electron_plotter->get1DHistogram("Top Layer " + std::to_string(layer) + " - y Residuals")->Fill(
                        y_residual); 
                } else {
                    positron_plotter->get2DHistogram("Top Layer " + std::to_string(layer) + " - Hit Positions")->Fill(
                        position[1], position[2]);
                    positron_plotter->get1DHistogram("Top Layer " + std::to_string(layer) + " - x Residuals")->Fill(
                        x_residual); 
                    positron_plotter->get1DHistogram("Top Layer " + std::to_string(layer) + " - y Residuals")->Fill(
                        y_residual); 
                }
            } else { 
                track_plotter->get2DHistogram("Bottom Layer " + std::to_string(layer) + " - Hit Positions")->Fill(
                        position[1], position[2]);
                track_plotter->get1DHistogram("Bottom Layer " + std::to_string(layer) + " - x Residuals")->Fill(
                        x_residual); 
                track_plotter->get1DHistogram("Bottom Layer " + std::to_string(layer) + " - y Residuals")->Fill(
                        y_residual); 
                if (track->getCharge() < 0) { 
                    electron_plotter->get2DHistogram("Bottom Layer " + std::to_string(layer) + " - Hit Positions")->Fill(
                        position[1], position[2]);
                    electron_plotter->get1DHistogram("Bottom Layer " + std::to_string(layer) + " - x Residuals")->Fill(
                        x_residual); 
                    electron_plotter->get1DHistogram("Bottom Layer " + std::to_string(layer) + " - y Residuals")->Fill(
                        y_residual); 
                } else {
                    positron_plotter->get2DHistogram("Bottom Layer " + std::to_string(layer) + " - Hit Positions")->Fill(
                        position[1], position[2]);
                    positron_plotter->get1DHistogram("Bottom Layer " + std::to_string(layer) + " - x Residuals")->Fill(
                        x_residual); 
                    positron_plotter->get1DHistogram("Bottom Layer " + std::to_string(layer) + " - y Residuals")->Fill(
                        y_residual); 
                }
            }
        }

        // FEE
        if (track->getCharge() < 0 && p_mag > 0.7) { 
            track_plotter->get1DHistogram("p - fee")->Fill(p_mag);
            track_plotter->get1DHistogram("px - fee")->Fill(p[0]);
            track_plotter->get1DHistogram("py - fee")->Fill(p[1]);
            track_plotter->get1DHistogram("pz - fee")->Fill(p[2]);
        }
    }

    if (event->getNumberOfTracks() == 2) {
     
        if (event->getTrack(0)->getCharge() < 0 && event->getTrack(1)->getCharge() > 0) {
            track_plotter->get2DHistogram("p[e+] v p[e-]")->Fill(
                    this->getMagnitude(event->getTrack(0)->getMomentum()),
                    this->getMagnitude(event->getTrack(1)->getMomentum()));
        } else if (event->getTrack(0)->getCharge() > 0 && event->getTrack(1)->getCharge() < 0) { 
            track_plotter->get2DHistogram("p[e+] v p[e-]")->Fill(
                    this->getMagnitude(event->getTrack(1)->getMomentum()),
                    this->getMagnitude(event->getTrack(0)->getMomentum()));
        } else if (event->getTrack(0)->getCharge() < 0 && event->getTrack(1)->getCharge() < 0) {
            track_plotter->get2DHistogram("p[e-] v p[e-]")->Fill(
                    this->getMagnitude(event->getTrack(0)->getMomentum()),
                    this->getMagnitude(event->getTrack(1)->getMomentum()));
        
            track_plotter->get2DHistogram("theta[e-] v theta[e-]")->Fill(
                    atan(event->getTrack(0)->getTanLambda()),
                    atan(event->getTrack(1)->getTanLambda()));
        } 
    }
}

void TrackAnalysis::finalize() {
    
    
    track_plotter->saveToRootFile("track_analysis.root");
    electron_plotter->saveToRootFile("electron_track_analysis.root");
    positron_plotter->saveToRootFile("positron_track_analysis.root");
    top_plotter->saveToRootFile("top_track_analysis.root");
    bottom_plotter->saveToRootFile("bottom_track_analysis.root");

    track_plotter->saveToPdf("track_analysis.pdf");
    electron_plotter->saveToPdf("electron_track_analysis.pdf");
    positron_plotter->saveToPdf("positron_track_analysis.pdf");
    top_plotter->saveToPdf("top_track_analysis.pdf");
    bottom_plotter->saveToPdf("bottom_track_analysis.pdf");
}

void TrackAnalysis::bookHistograms() { 

    // General event statistics
    track_plotter->setType("float"); 
    track_plotter->build1DHistogram("Number of tracks", 10, 0, 10);
    track_plotter->build1DHistogram("Track volume", 3, -1, 2);
    track_plotter->build1DHistogram("Hits per track", 6, 1, 7);
    track_plotter->build1DHistogram("Track charge", 3, -1, 2);
    track_plotter->build1DHistogram("chi2", 40, 0, 40);

    track_plotter->build1DHistogram("doca", 80, -10, 10);
    track_plotter->build1DHistogram("z0", 80, -2, 2);
    track_plotter->build1DHistogram("sin(phi0)", 40, -0.2, 0.2);
    track_plotter->build1DHistogram("curvature", 50, -0.001, 0.001);
    track_plotter->build1DHistogram("tan_lambda", 100, -0.1, 0.1);
    
    track_plotter->build1DHistogram("p", 50, 0, 2.0);
    track_plotter->build1DHistogram("pt", 50, -0.1, 0.2);
    track_plotter->build1DHistogram("px", 50, -0.1, 0.2); 
    track_plotter->build1DHistogram("py", 50, -0.15, 0.15); 
    track_plotter->build1DHistogram("pz", 50, 0, 2.0);
    track_plotter->build1DHistogram("p - fee", 50, .4, 2.0);
    track_plotter->build1DHistogram("px - fee", 50, -0.1, 0.2); 
    track_plotter->build1DHistogram("py - fee", 50, -0.15, 0.15); 
    track_plotter->build1DHistogram("pz - fee", 50, 0, 2.0);

    track_plotter->build2DHistogram("pz v px", 50, 0, 2.0, 50, -0.1, 0.2); 
    track_plotter->build2DHistogram("pz v py", 50, 0, 2.0, 50, -0.15, 0.15); 
    track_plotter->build2DHistogram("p v pt", 50, 0, 2.0, 50, -0.1, 0.2);

    track_plotter->build2DHistogram("track position at target", 80, -20, 20, 40, -5, 5);

    electron_plotter->setType("float")->setLineColor(kOrange + 9);
    electron_plotter->build1DHistogram("doca", 80, -10, 10);
    electron_plotter->build1DHistogram("z0", 80, -2, 2);
    electron_plotter->build1DHistogram("sin(phi0)", 40, -0.2, 0.2);
    electron_plotter->build1DHistogram("curvature", 50, -0.001, 0.001);
    electron_plotter->build1DHistogram("tan_lambda", 100, -0.1, 0.1);

    electron_plotter->build1DHistogram("p", 50, 0, 2.0);
    electron_plotter->build1DHistogram("pt", 50, -0.1, 0.2);
    electron_plotter->build1DHistogram("px", 50, -0.1, 0.2); 
    electron_plotter->build1DHistogram("py", 50, -0.15, 0.15); 
    electron_plotter->build1DHistogram("pz", 50, 0, 2.0); 
    electron_plotter->build1DHistogram("chi2", 40, 0, 40);
    electron_plotter->build1DHistogram("ep", 60, 0, 2);

    electron_plotter->build2DHistogram("track position at target", 80, -20, 20, 40, -5, 5); 

    positron_plotter->setType("float");
    positron_plotter->build1DHistogram("doca", 80, -10, 10);
    positron_plotter->build1DHistogram("z0", 80, -2, 2);
    positron_plotter->build1DHistogram("sin(phi0)", 40, -0.2, 0.2);
    positron_plotter->build1DHistogram("curvature", 50, -0.001, 0.001);
    positron_plotter->build1DHistogram("tan_lambda", 100, -0.1, 0.1);
            

    positron_plotter->build1DHistogram("p", 50, 0, 2.0);
    positron_plotter->build1DHistogram("pt", 50, -0.1, 0.2);
    positron_plotter->build1DHistogram("px", 50, -0.1, 0.2); 
    positron_plotter->build1DHistogram("py", 50, -0.15, 0.15); 
    positron_plotter->build1DHistogram("pz", 50, 0, 2.0); 
    positron_plotter->build1DHistogram("chi2", 40, 0, 40);

    positron_plotter->build2DHistogram("track position at target", 80, -20, 20, 40, -5, 5); 

    bottom_plotter->setType("float");
    top_plotter->build1DHistogram("doca", 80, -10, 10);
    top_plotter->build1DHistogram("z0", 80, -2, 2);
    top_plotter->build1DHistogram("sin(phi0)", 40, -0.2, 0.2);
    top_plotter->build1DHistogram("curvature", 50, -0.001, 0.001);
    top_plotter->build1DHistogram("tan_lambda", 100, -0.1, 0.1);

    top_plotter->build1DHistogram("p", 50, 0, 2.0);
    top_plotter->build1DHistogram("pt", 50, -0.1, 0.2);
    top_plotter->build1DHistogram("px", 50, -0.1, 0.2); 
    top_plotter->build1DHistogram("py", 50, -0.15, 0.15); 
    top_plotter->build1DHistogram("pz", 50, 0, 2.0); 
    top_plotter->build1DHistogram("chi2", 40, 0, 40);
    top_plotter->build1DHistogram("ep", 60, 0, 2);
    
    bottom_plotter->setType("float")->setLineColor(kOrange + 9);
    bottom_plotter->build1DHistogram("doca", 80, -10, 10);
    bottom_plotter->build1DHistogram("z0", 80, -2, 2);
    bottom_plotter->build1DHistogram("sin(phi0)", 40, -0.2, 0.2);
    bottom_plotter->build1DHistogram("curvature", 50, -0.001, 0.001);
    bottom_plotter->build1DHistogram("tan_lambda", 100, -0.1, 0.1);

    bottom_plotter->build1DHistogram("p", 50, 0, 2.0);
    bottom_plotter->build1DHistogram("pt", 50, -0.1, 0.2);
    bottom_plotter->build1DHistogram("px", 50, -0.1, 0.2); 
    bottom_plotter->build1DHistogram("py", 50, -0.15, 0.15); 
    bottom_plotter->build1DHistogram("pz", 50, 0, 2.0); 
    bottom_plotter->build1DHistogram("chi2", 40, 0, 40);
    bottom_plotter->build1DHistogram("ep", 60, 0, 2);
    
    for (int module_n = 1; module_n <= 6; ++module_n) { 
        track_plotter->build2DHistogram("Top Layer " + std::to_string(module_n) + " - Hit Positions",
                130, -100, 160, 40, -10, 70);
        track_plotter->build2DHistogram("Bottom Layer " + std::to_string(module_n) + " - Hit Positions",
                130, -100, 160, 40, -70, 10);
        electron_plotter->build2DHistogram("Top Layer " + std::to_string(module_n) + " - Hit Positions",
                130, -100, 160, 40, -10, 70);
        electron_plotter->build2DHistogram("Bottom Layer " + std::to_string(module_n) + " - Hit Positions",
                130, -100, 160, 40, -70, 10);
        positron_plotter->build2DHistogram("Top Layer " + std::to_string(module_n) + " - Hit Positions",
                130, -100, 160, 40, -10, 70);
        positron_plotter->build2DHistogram("Bottom Layer " + std::to_string(module_n) + " - Hit Positions",
                130, -100, 160, 40, -70, 10);
        track_plotter->build1DHistogram("Top Layer " + std::to_string(module_n) + " - x Residuals", 120, -60, 60); 
        track_plotter->build1DHistogram("Top Layer " + std::to_string(module_n) + " - y Residuals", 120, -60, 60); 
        track_plotter->build1DHistogram("Bottom Layer " + std::to_string(module_n) + " - x Residuals", 120, -60, 60); 
        track_plotter->build1DHistogram("Bottom Layer " + std::to_string(module_n) + " - y Residuals", 120, -60, 60); 
        electron_plotter->build1DHistogram("Top Layer " + std::to_string(module_n) + " - x Residuals", 120, -60, 60); 
        electron_plotter->build1DHistogram("Top Layer " + std::to_string(module_n) + " - y Residuals", 120, -60, 60); 
        electron_plotter->build1DHistogram("Bottom Layer " + std::to_string(module_n) + " - x Residuals", 120, -60, 60); 
        electron_plotter->build1DHistogram("Bottom Layer " + std::to_string(module_n) + " - y Residuals", 120, -60, 60); 
        positron_plotter->build1DHistogram("Top Layer " + std::to_string(module_n) + " - x Residuals", 120, -60, 60); 
        positron_plotter->build1DHistogram("Top Layer " + std::to_string(module_n) + " - y Residuals", 120, -60, 60); 
        positron_plotter->build1DHistogram("Bottom Layer " + std::to_string(module_n) + " - x Residuals", 120, -60, 60); 
        positron_plotter->build1DHistogram("Bottom Layer " + std::to_string(module_n) + " - y Residuals", 120, -60, 60); 
    }
    
    track_plotter->build2DHistogram("p[e+] v p[e-]", 50, 0, 2.0, 50, 0, 2.0);
    track_plotter->build2DHistogram("p[e-] v p[e-]", 50, 0, 2.0, 50, 0, 2.0);
    track_plotter->build2DHistogram("theta[e-] v theta[e-]", 100, -0.05, 0.05, 100, -0.05, 0.05);
}

std::string TrackAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;   
}

double TrackAnalysis::getMagnitude(std::vector<double> v) { 
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}
