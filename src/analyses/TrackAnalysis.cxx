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
      event_counter(0), 
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
    
    //if (!event->isPair1Trigger()) return;
    if (!event->isSingle1Trigger()) return;

    event_counter++; 

    track_plotter->get1DHistogram("Number of tracks")->Fill(event->getNumberOfTracks());

    if (event->getNumberOfTracks() == 0) { 
        std::cout << "Event " << event->getEventNumber() << " doesn't have any tracks." << std::endl;
    }

    if (event->getNumberOfTracks() == 1) { 
        track_plotter->getGraph("Track #chi^{2} vs Event")->SetPoint(event_counter - 1, 
                event->getEventNumber(), event->getTrack(0)->getChi2());

        track_plotter->getGraph("Track #Omega vs Event")->SetPoint(event_counter - 1, 
                event->getEventNumber(), event->getTrack(0)->getOmega());

        track_plotter->getGraph("Track sin(#phi_{0}) vs Event")->SetPoint(event_counter - 1, 
                event->getEventNumber(), sin(event->getTrack(0)->getPhi0()));
        
        track_plotter->getGraph("Track cos(#theta) vs Event")->SetPoint(event_counter - 1, 
                event->getEventNumber(), TrackExtrapolator::getCosTheta(event->getTrack(0))); 

        track_plotter->getGraph("Track d0 vs Event")->SetPoint(event_counter - 1, 
                event->getEventNumber(), event->getTrack(0)->getD0()); 
        
        track_plotter->getGraph("Track z0 vs Event")->SetPoint(event_counter - 1, 
                event->getEventNumber(), event->getTrack(0)->getZ0()); 
    }

    // Loop over all of the tracks in the event
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 

        // Get a track from the event
        track = event->getTrack(track_n); 

        //if (abs(track->getD0()) < 3) continue;
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
        track_plotter->get1DHistogram("cos(theta)")->Fill(TrackExtrapolator::getCosTheta(track));

        // Fill the track time parameters
        track_plotter->get1DHistogram("track time")->Fill(track->getTrackTime());
        track_plotter->get2DHistogram("track time v d0")->Fill(track->getTrackTime(), track->getD0());
        
        track_plotter->get2DHistogram("sin(phi0) v curvature")->Fill(sin(track->getPhi0()), track->getOmega()); 

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

        track_plotter->get2DHistogram("sin(phi0) v px")->Fill(sin(track->getPhi0()), p[0]); 
        track_plotter->get2DHistogram("sin(phi0) v py")->Fill(sin(track->getPhi0()), p[1]); 

        // Extrapolate the track to the target
        std::vector<double> track_pos_target = TrackExtrapolator::extrapolateTrack(track, 0);
        std::vector<double> track_pos_sp = TrackExtrapolator::extrapolateTrack(track, -100);
        
        track_plotter->get2DHistogram("track position at target")->Fill(track_pos_target[0], track_pos_target[1]);
        track_plotter->get2DHistogram("track position at sp")->Fill(track_pos_sp[0], track_pos_sp[1]);

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
            electron_plotter->get1DHistogram("cos(theta)")->Fill(TrackExtrapolator::getCosTheta(track));

            electron_plotter->get1DHistogram("track time")->Fill(track->getTrackTime());
            electron_plotter->get2DHistogram("track time v d0")->Fill(track->getTrackTime(), track->getD0());

            electron_plotter->get2DHistogram("sin(phi0) v curvature")->Fill(sin(track->getPhi0()), track->getOmega()); 
            electron_plotter->get2DHistogram("sin(phi0) v px")->Fill(sin(track->getPhi0()), p[0]); 
            electron_plotter->get2DHistogram("sin(phi0) v py")->Fill(sin(track->getPhi0()), p[1]); 

            if (track->getParticle()->getClusters()->GetEntriesFast() != 0) { 
               double energy = ((EcalCluster*) track->getParticle()->getClusters()->At(0))->getEnergy();
                electron_plotter->get1DHistogram("ep")->Fill(energy/p_mag);
            }

            electron_plotter->get2DHistogram("track position at target")->Fill(track_pos_target[0], track_pos_target[1]);
            electron_plotter->get2DHistogram("track position at sp")->Fill(track_pos_sp[0], track_pos_sp[1]);

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
            positron_plotter->get1DHistogram("cos(theta)")->Fill(TrackExtrapolator::getCosTheta(track));

            positron_plotter->get2DHistogram("sin(phi0) v curvature")->Fill(sin(track->getPhi0()), track->getOmega()); 
            positron_plotter->get2DHistogram("sin(phi0) v px")->Fill(sin(track->getPhi0()), p[0]); 
            positron_plotter->get2DHistogram("sin(phi0) v py")->Fill(sin(track->getPhi0()), p[1]); 

            positron_plotter->get1DHistogram("track time")->Fill(track->getTrackTime());
            positron_plotter->get2DHistogram("track time v d0")->Fill(track->getTrackTime(), track->getD0());

            positron_plotter->get2DHistogram("track position at target")->Fill(track_pos_target[0], track_pos_target[1]);
            positron_plotter->get2DHistogram("track position at sp")->Fill(track_pos_sp[0], track_pos_sp[1]);
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

            if (track->getCharge() < 0 && track->getParticle()->getClusters()->GetEntriesFast() != 0) { 
               double energy = ((EcalCluster*) track->getParticle()->getClusters()->At(0))->getEnergy();
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

            if (track->getCharge() < 0 && track->getParticle()->getClusters()->GetEntriesFast() != 0) { 
               double energy = ((EcalCluster*) track->getParticle()->getClusters()->At(0))->getEnergy();
                bottom_plotter->get1DHistogram("ep")->Fill(energy/p_mag);
            }
        }

        track_plotter->get1DHistogram("Hits per track")->Fill(track->getSvtHits()->GetEntriesFast()); 
    
        TRefArray* hits = track->getSvtHits();
        for (int hit_n = 0; hit_n < hits->GetSize(); ++hit_n) { 
            std::vector<double> position = ((SvtHit*) hits->At(hit_n))->getPosition();

            // Extrapolate the track to the target
            std::vector<double> track_pos_layer = TrackExtrapolator::extrapolateTrack(track, position[0]);

            double x_residual = track_pos_layer[1] - position[1]; 
            /*std::cout << "Track Pos: ( "
                      << track_pos_layer[0] << ", "
                      << track_pos_layer[1] << ", "
                      << track_pos_layer[2] << " )"
                      << std::endl;
            std::cout << "Hit Pos: ( " 
                      << position[0] << ", " 
                      << position[1] << ", " 
                      << position[2] << " )" 
                      << std::endl; */
            double y_residual = track_pos_layer[2] - position[2];

            int layer = ((SvtHit*) hits->At(hit_n))->getLayer();
            if (track->isTopTrack()) { 
                std::string top_name = "Top Layer " + std::to_string(layer);
                track_plotter->get2DHistogram(top_name + " - Hit Positions")->Fill(position[1], position[2]);
                track_plotter->get1DHistogram(top_name + " - x Residuals")->Fill(x_residual); 
                track_plotter->get1DHistogram(top_name + " - y Residuals")->Fill(y_residual); 
                if (track->getCharge() < 0) { 
                    electron_plotter->get2DHistogram(top_name + " - Hit Positions")->Fill(position[1], position[2]);
                    electron_plotter->get1DHistogram(top_name + " - x Residuals")->Fill(x_residual); 
                    electron_plotter->get1DHistogram(top_name + " - y Residuals")->Fill(y_residual); 
                } else {
                    positron_plotter->get2DHistogram(top_name + " - Hit Positions")->Fill(position[1], position[2]);
                    positron_plotter->get1DHistogram(top_name + " - x Residuals")->Fill(x_residual); 
                    positron_plotter->get1DHistogram(top_name + " - y Residuals")->Fill(y_residual); 
                }
            } else { 
                std::string bot_name = "Bottom Layer " + std::to_string(layer);
                track_plotter->get2DHistogram(bot_name + " - Hit Positions")->Fill(position[1], position[2]);
                track_plotter->get1DHistogram(bot_name + " - x Residuals")->Fill(x_residual); 
                track_plotter->get1DHistogram(bot_name + " - y Residuals")->Fill(y_residual); 
                if (track->getCharge() < 0) { 
                    electron_plotter->get2DHistogram(bot_name + " - Hit Positions")->Fill(position[1], position[2]);
                    electron_plotter->get1DHistogram(bot_name + " - x Residuals")->Fill(x_residual); 
                    electron_plotter->get1DHistogram(bot_name + " - y Residuals")->Fill(y_residual); 
                } else {
                    positron_plotter->get2DHistogram(bot_name + " - Hit Positions")->Fill(position[1], position[2]);
                    positron_plotter->get1DHistogram(bot_name + " - x Residuals")->Fill(x_residual); 
                    positron_plotter->get1DHistogram(bot_name + " - y Residuals")->Fill(y_residual); 
                }
            }

            track_plotter->get1DHistogram("track time - hit time")->Fill(track->getTrackTime() - ((SvtHit*) hits->At(hit_n))->getTime());
        
            if (track->getCharge() < 0) { 
                electron_plotter->get1DHistogram("track time - hit time")->Fill(track->getTrackTime() - ((SvtHit*) hits->At(hit_n))->getTime());
            } else { 
                positron_plotter->get1DHistogram("track time - hit time")->Fill(track->getTrackTime() - ((SvtHit*) hits->At(hit_n))->getTime());
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
    
        //if (abs(event->getTrack(0)->getD0()) > 3 || abs(event->getTrack(1)->getD0()) > 3) {

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
        //} 
    }
}

void TrackAnalysis::finalize() {
    

    for (int module_n = 1; module_n <= 6; ++module_n) {
        std::string top_name = "Top Layer " + std::to_string(module_n);
        std::string bot_name = "Bottom Layer " + std::to_string(module_n);
        
        // Residuals in x and y
        track_plotter->get1DHistogram(top_name + " - x Residuals")->Fit("gaus"); 
        track_plotter->get1DHistogram(top_name + " - y Residuals")->Fit("gaus"); 
        track_plotter->get1DHistogram(bot_name + " - x Residuals")->Fit("gaus"); 
        track_plotter->get1DHistogram(bot_name + " - y Residuals")->Fit("gaus"); 
        electron_plotter->get1DHistogram(top_name + " - x Residuals")->Fit("gaus"); 
        electron_plotter->get1DHistogram(top_name + " - y Residuals")->Fit("gaus"); 
        electron_plotter->get1DHistogram(bot_name + " - x Residuals")->Fit("gaus"); 
        electron_plotter->get1DHistogram(bot_name + " - y Residuals")->Fit("gaus"); 
        positron_plotter->get1DHistogram(top_name + " - x Residuals")->Fit("gaus"); 
        positron_plotter->get1DHistogram(top_name + " - y Residuals")->Fit("gaus"); 
        positron_plotter->get1DHistogram(bot_name + " - x Residuals")->Fit("gaus"); 
        positron_plotter->get1DHistogram(bot_name + " - y Residuals")->Fit("gaus"); 
    }

    
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
    track_plotter->build1DHistogram("Number of tracks", 10, 0, 10)->GetXaxis()->SetTitle("Number of tracks");
    track_plotter->build1DHistogram("Track volume", 3, -1, 2)->GetXaxis()->SetTitle("Track volume");
    track_plotter->build1DHistogram("Hits per track", 6, 1, 7)->GetXaxis()->SetTitle("Hits per track");
    track_plotter->build1DHistogram("Track charge", 3, -1, 2)->GetXaxis()->SetTitle("Track charge");
    track_plotter->build1DHistogram("chi2", 40, 0, 40)->GetXaxis()->SetTitle("Track #chi^{2}");

    track_plotter->build1DHistogram("doca", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    track_plotter->build1DHistogram("z0", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    track_plotter->build1DHistogram("sin(phi0)", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    track_plotter->build1DHistogram("curvature", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    track_plotter->build1DHistogram("tan_lambda", 100, -0.1, 0.1)->GetXaxis()->SetTitle("tan #lambda");
    track_plotter->build1DHistogram("cos(theta)", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");
    
    track_plotter->build1DHistogram("p", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    track_plotter->build1DHistogram("pt", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{t} [GeV]");
    track_plotter->build1DHistogram("px", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{x} [GeV]"); 
    track_plotter->build1DHistogram("py", 50, -0.15, 0.15)->GetXaxis()->SetTitle("p_{y} [GeV]"); 
    track_plotter->build1DHistogram("pz", 50, 0, 2.0)->GetXaxis()->SetTitle("p_{z} [GeV]");
    track_plotter->build1DHistogram("p - fee", 50, .4, 2.0)->GetXaxis()->SetTitle("FEE - p [GeV]");
    track_plotter->build1DHistogram("px - fee", 50, -0.1, 0.2)->GetXaxis()->SetTitle("FEE - p_{x} [GeV]"); 
    track_plotter->build1DHistogram("py - fee", 50, -0.15, 0.15)->GetXaxis()->SetTitle("FEE - p_{y} [GeV]"); 
    track_plotter->build1DHistogram("pz - fee", 50, 0, 2.0)->GetXaxis()->SetTitle("FEE - p_{z} [GeV]");

    track_plotter->build1DHistogram("track time", 100, -10, 10)->GetXaxis()->SetTitle("Track time [ns]");
    track_plotter->build1DHistogram("track time - hit time", 100, -10, 10)->GetXaxis()->SetTitle("Track time - Hit time [ns]");

    track_plotter->build2DHistogram("pz v px", 50, 0, 2.0, 50, -0.1, 0.2);
    track_plotter->get2DHistogram("pz v px")->GetXaxis()->SetTitle("p_{z} [GeV]"); 
    track_plotter->get2DHistogram("pz v px")->GetYaxis()->SetTitle("p_{x} [GeV]"); 
    track_plotter->build2DHistogram("pz v py", 50, 0, 2.0, 50, -0.15, 0.15); 
    track_plotter->get2DHistogram("pz v py")->GetXaxis()->SetTitle("p_{z} [GeV]"); 
    track_plotter->get2DHistogram("pz v py")->GetYaxis()->SetTitle("p_{y} [GeV]"); 
    track_plotter->build2DHistogram("p v pt", 50, 0, 2.0, 50, -0.1, 0.2);
    track_plotter->get2DHistogram("p v pt")->GetXaxis()->SetTitle("p [GeV]"); 
    track_plotter->get2DHistogram("p v pt")->GetYaxis()->SetTitle("p_{t} [GeV]"); 

    track_plotter->build2DHistogram("track position at target", 80, -20, 20, 40, -5, 5);
    track_plotter->get2DHistogram("track position at target")->GetXaxis()->SetTitle("x [mm]"); 
    track_plotter->get2DHistogram("track position at target")->GetYaxis()->SetTitle("y [mm]"); 
    track_plotter->build2DHistogram("track position at sp", 80, -20, 20, 80, -20, 20);
    track_plotter->get2DHistogram("track position at sp")->GetXaxis()->SetTitle("x [mm]"); 
    track_plotter->get2DHistogram("track position at sp")->GetYaxis()->SetTitle("y [mm]"); 
    track_plotter->build2DHistogram("sin(phi0) v px", 40, -0.2, 0.2, 50, -0.1, 0.2);
    track_plotter->get2DHistogram("sin(phi0) v px")->GetXaxis()->SetTitle("sin #phi_{0}"); 
    track_plotter->get2DHistogram("sin(phi0) v px")->GetYaxis()->SetTitle("p_{x} [GeV]"); 
    track_plotter->build2DHistogram("sin(phi0) v py", 40, -0.2, 0.2, 50, -0.15, 0.15);
    track_plotter->get2DHistogram("sin(phi0) v py")->GetXaxis()->SetTitle("sin #phi_{0}"); 
    track_plotter->get2DHistogram("sin(phi0) v py")->GetYaxis()->SetTitle("p_{y} [GeV]"); 
    track_plotter->build2DHistogram("sin(phi0) v curvature", 40, -0.2, 0.2, 50, -0.001, 0.001);
    track_plotter->get2DHistogram("sin(phi0) v curvature")->GetXaxis()->SetTitle("sin #phi_{0}"); 
    track_plotter->get2DHistogram("sin(phi0) v curvature")->GetYaxis()->SetTitle("#Omega"); 
    track_plotter->build2DHistogram("track time v d0", 100, -10, 10, 80, -10, 10);
    track_plotter->get2DHistogram("track time v d0")->GetXaxis()->SetTitle("Track time [ns]"); 
    track_plotter->get2DHistogram("track time v d0")->GetYaxis()->SetTitle("D0 [mm]"); 

    track_plotter->buildGraph("Track #chi^{2} vs Event")->GetXaxis()->SetTitle("Event Number");
    track_plotter->getGraph("Track #chi^{2} vs Event")->GetYaxis()->SetTitle("Track #chi^{2}");

    track_plotter->buildGraph("Track #Omega vs Event")->GetXaxis()->SetTitle("Event Number");
    track_plotter->getGraph("Track #Omega vs Event")->GetYaxis()->SetTitle("Track #Omega");

    track_plotter->buildGraph("Track sin(#phi_{0}) vs Event")->GetXaxis()->SetTitle("Event Number");
    track_plotter->getGraph("Track sin(#phi_{0}) vs Event")->GetYaxis()->SetTitle("Track sin(#phi_{0})");

    track_plotter->buildGraph("Track cos(#theta) vs Event")->GetXaxis()->SetTitle("Event Number");
    track_plotter->getGraph("Track cos(#theta) vs Event")->GetYaxis()->SetTitle("Track cos(#theta)");

    track_plotter->buildGraph("Track d0 vs Event")->GetXaxis()->SetTitle("Event Number");
    track_plotter->getGraph("Track d0 vs Event")->GetYaxis()->SetTitle("Track D0 (mm)");

    track_plotter->buildGraph("Track z0 vs Event")->GetXaxis()->SetTitle("Event Number");
    track_plotter->getGraph("Track z0 vs Event")->GetYaxis()->SetTitle("Track Z0 (mm)");


    electron_plotter->setType("float")->setLineColor(kOrange + 9);
    
    electron_plotter->build1DHistogram("doca", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    electron_plotter->build1DHistogram("z0", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    electron_plotter->build1DHistogram("sin(phi0)", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    electron_plotter->build1DHistogram("curvature", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    electron_plotter->build1DHistogram("tan_lambda", 100, -0.1, 0.1)->GetXaxis()->SetTitle("tan #lambda");
    electron_plotter->build1DHistogram("cos(theta)", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    electron_plotter->build1DHistogram("p", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    electron_plotter->build1DHistogram("pt", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{t} [GeV]");
    electron_plotter->build1DHistogram("px", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{x} [GeV]"); 
    electron_plotter->build1DHistogram("py", 50, -0.15, 0.15)->GetXaxis()->SetTitle("p_{y} [GeV]"); 
    electron_plotter->build1DHistogram("pz", 50, 0, 2.0)->GetXaxis()->SetTitle("p_{z} [GeV]");
    electron_plotter->build1DHistogram("chi2", 40, 0, 40)->GetXaxis()->SetTitle("Track #chi^{2}");
    electron_plotter->build1DHistogram("ep", 60, 0, 2)->GetXaxis()->SetTitle("E/p");

    electron_plotter->build1DHistogram("track time", 100, -10, 10)->GetXaxis()->SetTitle("Track time [ns]");
    electron_plotter->build1DHistogram("track time - hit time", 100, -10, 10)->GetXaxis()->SetTitle("Track time - Hit time [ns]");

    electron_plotter->build2DHistogram("track position at target", 80, -20, 20, 40, -5, 5);
    electron_plotter->get2DHistogram("track position at target")->GetXaxis()->SetTitle("x [mm]"); 
    electron_plotter->get2DHistogram("track position at target")->GetYaxis()->SetTitle("y [mm]"); 
    electron_plotter->build2DHistogram("track position at sp", 80, -20, 20, 80, -20, 20);
    electron_plotter->get2DHistogram("track position at sp")->GetXaxis()->SetTitle("x [mm]"); 
    electron_plotter->get2DHistogram("track position at sp")->GetYaxis()->SetTitle("y [mm]"); 
    electron_plotter->build2DHistogram("sin(phi0) v px", 40, -0.2, 0.2, 50, -0.1, 0.2);
    electron_plotter->get2DHistogram("sin(phi0) v px")->GetXaxis()->SetTitle("sin #phi_{0}"); 
    electron_plotter->get2DHistogram("sin(phi0) v px")->GetYaxis()->SetTitle("p_{x} [GeV]"); 
    electron_plotter->build2DHistogram("sin(phi0) v py", 40, -0.2, 0.2, 50, -0.15, 0.15);
    electron_plotter->get2DHistogram("sin(phi0) v py")->GetXaxis()->SetTitle("sin #phi_{0}"); 
    electron_plotter->get2DHistogram("sin(phi0) v py")->GetYaxis()->SetTitle("p_{y} [GeV]"); 
    electron_plotter->build2DHistogram("sin(phi0) v curvature", 40, -0.2, 0.2, 50, -0.001, 0.001);
    electron_plotter->get2DHistogram("sin(phi0) v curvature")->GetXaxis()->SetTitle("sin #phi_{0}"); 
    electron_plotter->get2DHistogram("sin(phi0) v curvature")->GetYaxis()->SetTitle("#Omega"); 
    electron_plotter->build2DHistogram("track time v d0", 100, -10, 10, 80, -10, 10);
    electron_plotter->get2DHistogram("track time v d0")->GetXaxis()->SetTitle("Track time [ns]"); 
    electron_plotter->get2DHistogram("track time v d0")->GetYaxis()->SetTitle("D0 [mm]"); 

    positron_plotter->build1DHistogram("doca", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    positron_plotter->build1DHistogram("z0", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    positron_plotter->build1DHistogram("sin(phi0)", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    positron_plotter->build1DHistogram("curvature", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    positron_plotter->build1DHistogram("tan_lambda", 100, -0.1, 0.1)->GetXaxis()->SetTitle("tan #lambda");
    positron_plotter->build1DHistogram("cos(theta)", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    positron_plotter->build1DHistogram("p", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    positron_plotter->build1DHistogram("pt", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{t} [GeV]");
    positron_plotter->build1DHistogram("px", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{x} [GeV]"); 
    positron_plotter->build1DHistogram("py", 50, -0.15, 0.15)->GetXaxis()->SetTitle("p_{y} [GeV]"); 
    positron_plotter->build1DHistogram("pz", 50, 0, 2.0)->GetXaxis()->SetTitle("p_{z} [GeV]");
    positron_plotter->build1DHistogram("chi2", 40, 0, 40)->GetXaxis()->SetTitle("Track #chi^{2}");

    positron_plotter->build1DHistogram("track time", 100, -10, 10)->GetXaxis()->SetTitle("Track time [ns]");
    positron_plotter->build1DHistogram("track time - hit time", 100, -10, 10)->GetXaxis()->SetTitle("Track time - Hit time [ns]");

    positron_plotter->build2DHistogram("track position at target", 80, -20, 20, 40, -5, 5);
    positron_plotter->get2DHistogram("track position at target")->GetXaxis()->SetTitle("x [mm]"); 
    positron_plotter->get2DHistogram("track position at target")->GetYaxis()->SetTitle("y [mm]"); 
    positron_plotter->build2DHistogram("track position at sp", 80, -20, 20, 80, -20, 20);
    positron_plotter->get2DHistogram("track position at sp")->GetXaxis()->SetTitle("x [mm]"); 
    positron_plotter->get2DHistogram("track position at sp")->GetYaxis()->SetTitle("y [mm]"); 
    positron_plotter->build2DHistogram("sin(phi0) v px", 40, -0.2, 0.2, 50, -0.1, 0.2);
    positron_plotter->get2DHistogram("sin(phi0) v px")->GetXaxis()->SetTitle("sin #phi_{0}"); 
    positron_plotter->get2DHistogram("sin(phi0) v px")->GetYaxis()->SetTitle("p_{x} [GeV]"); 
    positron_plotter->build2DHistogram("sin(phi0) v py", 40, -0.2, 0.2, 50, -0.15, 0.15);
    positron_plotter->get2DHistogram("sin(phi0) v py")->GetXaxis()->SetTitle("sin #phi_{0}"); 
    positron_plotter->get2DHistogram("sin(phi0) v py")->GetYaxis()->SetTitle("p_{y} [GeV]"); 
    positron_plotter->build2DHistogram("sin(phi0) v curvature", 40, -0.2, 0.2, 50, -0.001, 0.001);
    positron_plotter->get2DHistogram("sin(phi0) v curvature")->GetXaxis()->SetTitle("sin #phi_{0}"); 
    positron_plotter->get2DHistogram("sin(phi0) v curvature")->GetYaxis()->SetTitle("#Omega"); 
    positron_plotter->build2DHistogram("track time v d0", 100, -10, 10, 80, -10, 10);
    positron_plotter->get2DHistogram("track time v d0")->GetXaxis()->SetTitle("Track time [ns]"); 
    positron_plotter->get2DHistogram("track time v d0")->GetYaxis()->SetTitle("D0 [mm]"); 

    top_plotter->setType("float");
    top_plotter->build1DHistogram("doca", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    top_plotter->build1DHistogram("z0", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    top_plotter->build1DHistogram("sin(phi0)", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    top_plotter->build1DHistogram("curvature", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    top_plotter->build1DHistogram("tan_lambda", 100, -0.1, 0.1)->GetXaxis()->SetTitle("tan #lambda");
    top_plotter->build1DHistogram("cos(theta)", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    top_plotter->build1DHistogram("p", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    top_plotter->build1DHistogram("pt", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{t} [GeV]");
    top_plotter->build1DHistogram("px", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{x} [GeV]"); 
    top_plotter->build1DHistogram("py", 50, -0.15, 0.15)->GetXaxis()->SetTitle("p_{y} [GeV]"); 
    top_plotter->build1DHistogram("pz", 50, 0, 2.0)->GetXaxis()->SetTitle("p_{z} [GeV]");
    top_plotter->build1DHistogram("chi2", 40, 0, 40)->GetXaxis()->SetTitle("Track #chi^{2}");
    top_plotter->build1DHistogram("ep", 60, 0, 2)->GetXaxis()->SetTitle("E/p");

    bottom_plotter->build1DHistogram("doca", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    bottom_plotter->build1DHistogram("z0", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    bottom_plotter->build1DHistogram("sin(phi0)", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    bottom_plotter->build1DHistogram("curvature", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    bottom_plotter->build1DHistogram("tan_lambda", 100, -0.1, 0.1)->GetXaxis()->SetTitle("tan #lambda");
    bottom_plotter->build1DHistogram("cos(theta)", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    bottom_plotter->build1DHistogram("p", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    bottom_plotter->build1DHistogram("pt", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{t} [GeV]");
    bottom_plotter->build1DHistogram("px", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{x} [GeV]"); 
    bottom_plotter->build1DHistogram("py", 50, -0.15, 0.15)->GetXaxis()->SetTitle("p_{y} [GeV]"); 
    bottom_plotter->build1DHistogram("pz", 50, 0, 2.0)->GetXaxis()->SetTitle("p_{z} [GeV]");
    bottom_plotter->build1DHistogram("chi2", 40, 0, 40)->GetXaxis()->SetTitle("Track #chi^{2}");
    bottom_plotter->build1DHistogram("ep", 60, 0, 2)->GetXaxis()->SetTitle("E/p");

    double range = 0.5;    
    for (int module_n = 1; module_n <= 6; ++module_n) {
        std::string top_name = "Top Layer " + std::to_string(module_n);
        std::string bot_name = "Bottom Layer " + std::to_string(module_n);

        // Track positions at the target
        track_plotter->build2DHistogram(top_name + " - Hit Positions", 130, -100, 160, 40, -10, 70);
        track_plotter->get2DHistogram(top_name + " - Hit Positions")->GetXaxis()->SetTitle("Top Hit x [mm]"); 
        track_plotter->get2DHistogram(top_name + " - Hit Positions")->GetYaxis()->SetTitle("Top Hit y [mm]"); 
        track_plotter->build2DHistogram(bot_name + " - Hit Positions", 130, -100, 160, 40, -70, 10);
        track_plotter->get2DHistogram(bot_name + " - Hit Positions")->GetXaxis()->SetTitle("Bottom Hit x [mm]"); 
        track_plotter->get2DHistogram(bot_name + " - Hit Positions")->GetYaxis()->SetTitle("Bottom Hit y [mm]");
        electron_plotter->build2DHistogram(top_name + " - Hit Positions", 130, -100, 160, 40, -10, 70);
        electron_plotter->get2DHistogram(top_name + " - Hit Positions")->GetXaxis()->SetTitle("Top Hit x [mm]"); 
        electron_plotter->get2DHistogram(top_name + " - Hit Positions")->GetYaxis()->SetTitle("Top Hit y [mm]"); 
        electron_plotter->build2DHistogram(bot_name + " - Hit Positions", 130, -100, 160, 40, -70, 10);
        electron_plotter->get2DHistogram(bot_name + " - Hit Positions")->GetXaxis()->SetTitle("Bottom Hit x [mm]"); 
        electron_plotter->get2DHistogram(bot_name + " - Hit Positions")->GetYaxis()->SetTitle("Bottom Hit y [mm]"); 
        positron_plotter->build2DHistogram(top_name + " - Hit Positions", 130, -100, 160, 40, -10, 70);
        positron_plotter->get2DHistogram(top_name + " - Hit Positions")->GetXaxis()->SetTitle("Top Hit x [mm]"); 
        positron_plotter->get2DHistogram(top_name + " - Hit Positions")->GetYaxis()->SetTitle("Top Hit y [mm]"); 
        positron_plotter->build2DHistogram(bot_name + " - Hit Positions", 130, -100, 160, 40, -70, 10);
        positron_plotter->get2DHistogram(bot_name + " - Hit Positions")->GetXaxis()->SetTitle("Bottom Hit x [mm]"); 
        positron_plotter->get2DHistogram(bot_name + " - Hit Positions")->GetYaxis()->SetTitle("Bottom Hit y [mm]"); 

        // Residuals in x and y
        
        track_plotter->build1DHistogram(top_name + " - x Residuals", 200, -module_n*range, module_n*range)->GetXaxis()->SetTitle("Top Layer x Residual"); 
        track_plotter->build1DHistogram(top_name + " - y Residuals", 200, -module_n*range, module_n*range)->GetXaxis()->SetTitle("Top Layer y Residual"); 
        track_plotter->build1DHistogram(bot_name + " - x Residuals", 200, -module_n*range, module_n*range)->GetXaxis()->SetTitle("Bottom Layer x Residual"); 
        track_plotter->build1DHistogram(bot_name + " - y Residuals", 200, -module_n*range, module_n*range)->GetXaxis()->SetTitle("Bottom Layer y Residual"); 
        track_plotter->build1DHistogram(top_name + " - Number of hits per cluster", 5, 0, 5); 
        track_plotter->build1DHistogram(bot_name + " - Number of hits per cluster", 5, 0, 5); 
        electron_plotter->build1DHistogram(top_name + " - x Residuals", 200, -module_n*range, module_n*range)->GetXaxis()->SetTitle("Top Layer x Residual"); 
        electron_plotter->build1DHistogram(top_name + " - y Residuals", 200, -module_n*range, module_n*range)->GetXaxis()->SetTitle("Top Layer y Residual"); 
        electron_plotter->build1DHistogram(bot_name + " - x Residuals", 200, -module_n*range, module_n*range)->GetXaxis()->SetTitle("Bottom Layer x Residual"); 
        electron_plotter->build1DHistogram(bot_name + " - y Residuals", 200, -module_n*range, module_n*range)->GetXaxis()->SetTitle("Bottom Layer y Residual"); 
        positron_plotter->build1DHistogram(top_name + " - x Residuals", 200, -module_n*range, module_n*range)->GetXaxis()->SetTitle("Top Layer x Residual"); 
        positron_plotter->build1DHistogram(top_name + " - y Residuals", 200, -module_n*range, module_n*range)->GetXaxis()->SetTitle("Top Layer y Residual"); 
        positron_plotter->build1DHistogram(bot_name + " - x Residuals", 200, -module_n*range, module_n*range)->GetXaxis()->SetTitle("Bottom Layer x Residual"); 
        positron_plotter->build1DHistogram(bot_name + " - y Residuals", 200, -module_n*range, module_n*range)->GetXaxis()->SetTitle("Bottom Layer y Residual"); 
    }
    
    track_plotter->build2DHistogram("p[e+] v p[e-]", 50, 0, 2.0, 50, 0, 2.0);
    track_plotter->get2DHistogram("p[e+] v p[e-]")->GetXaxis()->SetTitle("p[e+] [GeV]"); 
    track_plotter->get2DHistogram("p[e+] v p[e-]")->GetYaxis()->SetTitle("p[e-] [GeV]"); 
    track_plotter->build2DHistogram("p[e-] v p[e-]", 50, 0, 2.0, 50, 0, 2.0);
    track_plotter->get2DHistogram("p[e-] v p[e-]")->GetXaxis()->SetTitle("p[e-] [GeV]"); 
    track_plotter->get2DHistogram("p[e-] v p[e-]")->GetYaxis()->SetTitle("p[e-] [GeV]"); 
    track_plotter->build2DHistogram("theta[e-] v theta[e-]", 100, -0.05, 0.05, 100, -0.05, 0.05);
    track_plotter->get2DHistogram("theta[e-] v theta[e-]")->GetXaxis()->SetTitle("theta[e-]"); 
    track_plotter->get2DHistogram("theta[e-] v theta[e-]")->GetYaxis()->SetTitle("theta[e-]"); 
}

std::string TrackAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;   
}

double TrackAnalysis::getMagnitude(std::vector<double> v) { 
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}
