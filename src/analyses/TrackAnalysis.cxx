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
      canvas(NULL),
      output_file(new TFile("track_analysis_results.root", "RECREATE")),
      electron_plotter(new Plotter()),
      class_name("TrackAnalysis") {  
}

TrackAnalysis::~TrackAnalysis() { 
    delete canvas;
    delete output_file;
    delete electron_plotter;
}

void TrackAnalysis::initialize() { 

    this->bookHistograms();
}

void TrackAnalysis::processEvent(HpsEvent* event) {

    if (!event->isPair1Trigger()) return;

    track_plots["Number of tracks"]->Fill(event->getNumberOfTracks());

    // Loop over all of the tracks in the event
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 

        // Get a track from the event
        track = event->getTrack(track_n); 

        // Fill the plots providing general event information related to tracks
        int track_volume = track->isTopTrack() ? 0 : 1; 
        track_plots["Track volume"]->Fill(track_volume);
        track_plots["Track charge"]->Fill(track->getCharge());
        track_plots["chi2"]->Fill(track->getChi2());

        // Fill the track parameter plots for all tracks
        track_plots["doca"]->Fill(track->getD0());
        track_plots["z0"]->Fill(track->getZ0());
        track_plots["sin(phi0)"]->Fill(sin(track->getPhi0()));
        track_plots["curvature"]->Fill(track->getOmega());
        track_plots["tan_lambda"]->Fill(track->getTanLambda());
       
         
        // Calculate the momentum magnitude and transverse momentum
        std::vector<double> p = track->getMomentum();
        double p_mag = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        double pt = sqrt(p[0]*p[0] + p[1]*p[1]);
        
        // Fill the momentum plots
        track_plots["p"]->Fill(p_mag);
        track_plots["pt"]->Fill(pt);
        track_plots["px"]->Fill(p[0]);
        track_plots["py"]->Fill(p[1]);
        track_plots["pz"]->Fill(p[2]);
        track_plots_2d["pz v px"]->Fill(p[2], p[0]);
        track_plots_2d["pz v py"]->Fill(p[2], p[1]);
        track_plots_2d["p v pt"]->Fill(p_mag, pt);

        // Fill the electron and positron plots
        if (track->getCharge() < 0) { 
            electron_track_plots["p"]->Fill(p_mag);
            electron_track_plots["pt"]->Fill(pt);
            electron_track_plots["px"]->Fill(p[0]);
            electron_track_plots["py"]->Fill(p[1]);
            electron_track_plots["pz"]->Fill(p[2]);
            electron_track_plots["chi2"]->Fill(track->getChi2());
            electron_track_plots["doca"]->Fill(track->getD0());

            electron_plotter->get1DHistogram("doca")->Fill(track->getD0());

            electron_track_plots["z0"]->Fill(track->getZ0());
            electron_track_plots["sin(phi0)"]->Fill(sin(track->getPhi0()));
            electron_track_plots["curvature"]->Fill(track->getOmega());
            electron_track_plots["tan_lambda"]->Fill(track->getTanLambda());

            if (((HpsParticle*) track->getParticle().GetObject())->getClusters()->GetEntriesFast() != 0) { 
               double energy = ((EcalCluster*) ((HpsParticle*) track->getParticle().GetObject())->getClusters()->At(0))->getEnergy();
                electron_track_plots["ep"]->Fill(energy/p_mag);
            }
        } else { 
            positron_track_plots["p"]->Fill(p_mag);
            positron_track_plots["pt"]->Fill(pt);
            positron_track_plots["px"]->Fill(p[0]);
            positron_track_plots["py"]->Fill(p[1]);
            positron_track_plots["pz"]->Fill(p[2]);
            positron_track_plots["chi2"]->Fill(track->getChi2());
            positron_track_plots["doca"]->Fill(track->getD0());
            positron_track_plots["z0"]->Fill(track->getZ0());
            positron_track_plots["sin(phi0)"]->Fill(sin(track->getPhi0()));
            positron_track_plots["curvature"]->Fill(track->getOmega());
            positron_track_plots["tan_lambda"]->Fill(track->getTanLambda());
        }
    
        // Fill the top and bottom track plots
        if (track->isTopTrack()) { 
            
            top_track_plots["p"]->Fill(p_mag);
            top_track_plots["pt"]->Fill(pt);
            top_track_plots["px"]->Fill(p[0]);
            top_track_plots["py"]->Fill(p[1]);
            top_track_plots["pz"]->Fill(p[2]);
            top_track_plots["chi2"]->Fill(track->getChi2());
            top_track_plots["doca"]->Fill(track->getD0());
            top_track_plots["z0"]->Fill(track->getZ0());
            top_track_plots["sin(phi0)"]->Fill(sin(track->getPhi0()));
            top_track_plots["curvature"]->Fill(track->getOmega());
            top_track_plots["tan_lambda"]->Fill(track->getTanLambda());

            if (track->getCharge() < 0 && ((HpsParticle*) track->getParticle().GetObject())->getClusters()->GetEntriesFast() != 0) { 
               double energy = ((EcalCluster*) ((HpsParticle*) track->getParticle().GetObject())->getClusters()->At(0))->getEnergy();
                top_track_plots["ep"]->Fill(energy/p_mag);
            }
        
        } else { 
            
            bot_track_plots["p"]->Fill(p_mag);
            bot_track_plots["pt"]->Fill(pt);
            bot_track_plots["px"]->Fill(p[0]);
            bot_track_plots["py"]->Fill(p[1]);
            bot_track_plots["pz"]->Fill(p[2]);
            bot_track_plots["chi2"]->Fill(track->getChi2());
            bot_track_plots["doca"]->Fill(track->getD0());
            bot_track_plots["z0"]->Fill(track->getZ0());
            bot_track_plots["sin(phi0)"]->Fill(sin(track->getPhi0()));
            bot_track_plots["curvature"]->Fill(track->getOmega());
            bot_track_plots["tan_lambda"]->Fill(track->getTanLambda());

            if (track->getCharge() < 0 && ((HpsParticle*) track->getParticle().GetObject())->getClusters()->GetEntriesFast() != 0) { 
               double energy = ((EcalCluster*) ((HpsParticle*) track->getParticle().GetObject())->getClusters()->At(0))->getEnergy();
                bot_track_plots["ep"]->Fill(energy/p_mag);
            }
        }

        track_plots["Hits per track"]->Fill(track->getSvtHits()->GetEntriesFast()); 
    
        TRefArray* hits = track->getSvtHits();
        for (int hit_n = 0; hit_n < hits->GetSize(); ++hit_n) { 
            std::vector<double> position = ((SvtHit*) hits->At(hit_n))->getPosition();
            int layer = ((SvtHit*) hits->At(hit_n))->getLayer();
            hit_position_plots["Top Module " + std::to_string(layer)]->Fill(position[1], position[2]);
        }

        // FEE
        if (track->getCharge() < 0 && p_mag > 0.7) { 
            track_plots["p - fee"]->Fill(p_mag);
            track_plots["px - fee"]->Fill(p[0]);
            track_plots["py - fee"]->Fill(p[1]);
            track_plots["pz - fee"]->Fill(p[2]);
        }
    }

    if (event->getNumberOfTracks() == 2) { 
     
        if (event->getTrack(0)->getCharge() < 0 && event->getTrack(1)->getCharge() > 0) {
            track_epem_2d_plots["p[e+] v p[e-]"]->Fill(
                    this->getMagnitude(event->getTrack(0)->getMomentum()),
                    this->getMagnitude(event->getTrack(1)->getMomentum()));
        } else if (event->getTrack(0)->getCharge() > 0 && event->getTrack(1)->getCharge() < 0) { 
            track_epem_2d_plots["p[e+] v p[e-]"]->Fill(
                    this->getMagnitude(event->getTrack(1)->getMomentum()),
                    this->getMagnitude(event->getTrack(0)->getMomentum()));
        } else if (event->getTrack(0)->getCharge() < 0 && event->getTrack(1)->getCharge() < 0) {
            track_emem_2d_plots["p[e-] v p[e-]"]->Fill(
                    this->getMagnitude(event->getTrack(0)->getMomentum()),
                    this->getMagnitude(event->getTrack(1)->getMomentum()));
        
            track_emem_2d_plots["theta[e-] v theta[e-]"]->Fill(
                    atan(event->getTrack(0)->getTanLambda()),
                    atan(event->getTrack(1)->getTanLambda()));
        } 
    }
}

void TrackAnalysis::finalize() {

    track_plots["p - fee"]->Fit("gaus");

    canvas->Print("track_analysis_results.pdf[");

    std::unordered_map<std::string, TH1F*>::iterator track_plots_it = track_plots.begin();
    for (track_plots_it; track_plots_it != track_plots.end(); ++track_plots_it) { 
        track_plots_it->second->Draw();
        track_plots_it->second->Write();
        canvas->Print("track_analysis_results.pdf(");
    } 

    std::unordered_map<std::string, TH1F*>::iterator electron_track_plots_it = electron_track_plots.begin();
    for (electron_track_plots_it; 
            electron_track_plots_it != electron_track_plots.end();
            ++electron_track_plots_it) { 
        electron_track_plots_it->second->Draw();
        electron_track_plots_it->second->Write();
        
        if (positron_track_plots[electron_track_plots_it->first] != NULL) { 
            positron_track_plots[electron_track_plots_it->first]->Draw("same");
            positron_track_plots[electron_track_plots_it->first]->Write();
        }
        canvas->Print("track_analysis_results.pdf(");
    }

    std::unordered_map<std::string, TH1F*>::iterator bot_track_plots_it = bot_track_plots.begin();
    for (bot_track_plots_it; 
            bot_track_plots_it != bot_track_plots.end();
            ++bot_track_plots_it) { 
        bot_track_plots_it->second->Draw();
        bot_track_plots_it->second->Write();
        
        if (top_track_plots[bot_track_plots_it->first] != NULL) { 
            top_track_plots[bot_track_plots_it->first]->Draw("same");
            top_track_plots[bot_track_plots_it->first]->Write();
        }
        canvas->Print("track_analysis_results.pdf(");
    }

    std::unordered_map<std::string, TH2F*>::iterator hit_position_plots_it = hit_position_plots.begin();
    for (hit_position_plots_it;
            hit_position_plots_it != hit_position_plots.end(); 
            ++hit_position_plots_it) { 
        hit_position_plots_it->second->Draw("colz");
        canvas->Print("track_analysis_results.pdf("); 
    }
   
    std::unordered_map<std::string, TH2F*>::iterator track_plots_2d_it = track_plots_2d.begin();
    for (track_plots_2d_it;
            track_plots_2d_it != track_plots_2d.end(); 
            ++track_plots_2d_it) { 
        track_plots_2d_it->second->Draw("colz");
        canvas->Print("track_analysis_results.pdf("); 
    }

    track_epem_2d_plots["p[e+] v p[e-]"]->Draw("colz");
    canvas->Print("track_analysis_results.pdf("); 

    track_emem_2d_plots["p[e-] v p[e-]"]->Draw("colz");
    canvas->Print("track_analysis_results.pdf("); 

    track_emem_2d_plots["theta[e-] v theta[e-]"]->Draw("colz");
    canvas->Print("track_analysis_results.pdf("); 
    
    electron_plotter->get1DHistogram("doca")->Draw();
    canvas->Print("track_analysis_results.pdf("); 

    canvas->Print("track_analysis_results.pdf]");
    output_file->Close(); 
}

void TrackAnalysis::bookHistograms() { 

    canvas = new TCanvas("canvas", "canvas", 500, 500);
  
    // General event statistics 
    track_plots["Number of tracks"] = new TH1F("number_of_tracks", "number_of_tracks", 10, 0, 10);
    track_plots["Track volume"] = new TH1F("track_volume", "track_volume", 3, -1, 2);
    track_plots["Hits per track"] = new TH1F("hits_per_track", "hits_per_track", 6, 1, 7);
    track_plots["Track charge"] = new TH1F("track_charge", "track_charge", 3, -1, 2);
    track_plots["chi2"] = new TH1F("chi2", "chi2", 40, 0, 40);

    track_plots["doca"] = new TH1F("doca", "doca", 40, -20, 20);
    track_plots["z0"] = new TH1F("z0", "z0", 100, -5, 5);
    track_plots["sin(phi0)"] = new TH1F("sin(phi0)", "sin(phi0)", 40, -0.4, 0.4);
    track_plots["curvature"] = new TH1F("curvature", "curvature", 25, -0.0025, 0.0025);
    track_plots["tan_lambda"] = new TH1F("tan_lambda", "tan_lambda", 100, -0.1, 0.1);
    
    track_plots["p"] = new TH1F("p", "p", 50, 0, 2.0);
    track_plots["pt"] = new TH1F("pt", "pt", 50, -0.1, 0.2);
    track_plots["px"] = new TH1F("px", "px", 50, -0.1, 0.2); 
    track_plots["py"] = new TH1F("py", "py", 50, -0.15, 0.15); 
    track_plots["pz"] = new TH1F("pz", "pz", 50, 0, 2.0);
    track_plots["p - fee"] = new TH1F("p_fee", "p_fee", 50, .4, 2.0);
    track_plots["px - fee"] = new TH1F("px_fee", "px_fee", 50, -0.1, 0.2); 
    track_plots["py - fee"] = new TH1F("py_fee", "py_fee", 50, -0.15, 0.15); 
    track_plots["pz - fee"] = new TH1F("pz_fee", "pz_fee", 50, 0, 2.0);

    track_plots_2d["pz v px"] = new TH2F("pz_v_px", "pz_v_px", 50, 0, 2.0, 50, -0.1, 0.2); 
    track_plots_2d["pz v py"] = new TH2F("pz_v_py", "pz_v_py", 50, 0, 2.0, 50, -0.1, 0.2); 
    track_plots_2d["p v pt"] = new TH2F("p_v_pt", "p_v_pt", 50, 0, 2.0, 50, -0.1, 0.2); 

    electron_track_plots["doca"] = new TH1F("doca_electron", "doca_electron", 40, -20, 20);
    electron_track_plots["z0"] = new TH1F("z0_electron", "z0_electron", 100, -5, 5);
    electron_track_plots["sin(phi0)"] = new TH1F("sin(phi0)_electron", "sin(phi0)_electron", 40, -0.4, 0.4);
    electron_track_plots["curvature"] = new TH1F("curvature_electron", "curvature_electron", 25, -0.0025, 0.0025);
    electron_track_plots["tan_lambda"] = new TH1F("tan_lambda_electron", "tan_lambda_electron", 100, -0.1, 0.1);

    electron_plotter->setType("float")->setLineColor(kRed)->build1DHistogram("doca", 40, -20, 20);

    electron_track_plots["p"] = new TH1F("p_electron", "p_electron", 50, 0, 2.0);
    electron_track_plots["pt"] = new TH1F("pt_electron", "pt_electron", 50, -0.1, 0.2);
    electron_track_plots["px"] = new TH1F("px_electron", "px_electron", 50, -0.1, 0.2); 
    electron_track_plots["py"] = new TH1F("py_electron", "py_electron", 50, -0.2, 0.2); 
    electron_track_plots["pz"] = new TH1F("pz_electron", "pz_electron", 50, 0, 3.0); 
    electron_track_plots["chi2"] = new TH1F("chi2_electron", "chi2_electron", 40, 0, 40);
    electron_track_plots["ep"] = new TH1F("ep", "ep", 60, 0, 2);

    positron_track_plots["doca"] = new TH1F("doca_positron", "doca_positron", 40, -20, 20);
    positron_track_plots["z0"] = new TH1F("z0_positron", "z0_positron", 100, -5, 5);
    positron_track_plots["sin(phi0)"] = new TH1F("sin(phi0)_positron", "sin(phi0)_positron", 60, -60, 60);
    positron_track_plots["curvature"] = new TH1F("curvature_positron", "curvature_positron", 25, -0.0025, 0.0025);
    positron_track_plots["tan_lambda"] = new TH1F("tan_lambda_positron", "tan_lambda_positron", 100, -0.1, 0.1);

    positron_track_plots["p"] = new TH1F("p_positron", "p_positron", 50, 0, 2.0);
    positron_track_plots["pt"] = new TH1F("pt_positron", "pt_positron", 50, -0.1, 0.2);
    positron_track_plots["px"] = new TH1F("px_positron", "px_positron", 50, -0.1, 0.2); 
    positron_track_plots["py"] = new TH1F("py_positron", "py_positron", 50, -0.15, 0.15); 
    positron_track_plots["pz"] = new TH1F("pz_positron", "pz_positron", 50, 0, 3.0); 
    positron_track_plots["chi2"] = new TH1F("chi2_positron", "chi2_positron", 40, 0, 40);

    top_track_plots["doca"] = new TH1F("doca_top", "doca_top", 40, -20, 20);
    top_track_plots["z0"] = new TH1F("z0_top", "z0_top", 100, -5, 5);
    top_track_plots["sin(phi0)"] = new TH1F("sin(phi0)_top", "sin(phi0)_top", 60, -60, 60);
    top_track_plots["curvature"] = new TH1F("curvature_top", "curvature_top", 25, -0.0025, 0.0025);
    top_track_plots["tan_lambda"] = new TH1F("tan_lambda_top", "tan_lambda_top", 100, -0.1, 0.1);

    top_track_plots["p"] = new TH1F("p_top", "p_top", 50, 0, 2.0);
    top_track_plots["pt"] = new TH1F("pt_top", "pt_top", 50, -0.1, 0.2);
    top_track_plots["px"] = new TH1F("px_top", "px_top", 50, -0.1, 0.2); 
    top_track_plots["py"] = new TH1F("py_top", "py_top", 50, -0.2, 0.2); 
    top_track_plots["pz"] = new TH1F("pz_top", "pz_top", 50, 0, 3.0); 
    top_track_plots["chi2"] = new TH1F("chi2_top", "chi2_top", 40, 0, 40);
    top_track_plots["ep"] = new TH1F("ep_top", "ep_top", 60, 0, 2);

    bot_track_plots["doca"] = new TH1F("doca_bot", "doca_bot", 40, -20, 20);
    bot_track_plots["z0"] = new TH1F("z0_bot", "z0_bot", 100, -5, 5);
    bot_track_plots["sin(phi0)"] = new TH1F("sin(phi0)_bot", "sin(phi0)_bot", 60, -60, 60);
    bot_track_plots["curvature"] = new TH1F("curvature_bot", "curvature_bot", 25, -0.0025, 0.0025);
    bot_track_plots["tan_lambda"] = new TH1F("tan_lambda_bot", "tan_lambda_bot", 100, -0.1, 0.1);

    bot_track_plots["p"] = new TH1F("p_bot", "p_bot", 50, 0, 2.0);
    bot_track_plots["pt"] = new TH1F("pt_bot", "pt_bot", 50, -0.1, 0.2);
    bot_track_plots["px"] = new TH1F("px_bot", "px_bot", 50, -0.1, 0.2); 
    bot_track_plots["py"] = new TH1F("py_bot", "py_bot", 50, -0.2, 0.2); 
    bot_track_plots["pz"] = new TH1F("pz_bot", "pz_bot", 50, 0, 3.0); 
    bot_track_plots["chi2"] = new TH1F("chi2_bot", "chi2_bot", 40, 0, 40);
    bot_track_plots["ep"] = new TH1F("ep_bot", "ep_bot", 60, 0, 2);


    for (int module_n = 1; module_n <= 6; ++module_n) { 
        hit_position_plots["Top Module " + std::to_string(module_n)] 
            = new TH2F(("top_module_" + std::to_string(module_n)).c_str(),
                       ("top_module_" + std::to_string(module_n)).c_str(), 
                       50, -100, 100, 50, 0, 50);
        hit_position_plots["Bottom Module " + std::to_string(module_n)] 
            = new TH2F(("bottom_module_" + std::to_string(module_n)).c_str(),
                       ("bottom_module_" + std::to_string(module_n)).c_str(), 
                       50, -100, 100, 50, 0, 50); 
    }

    track_epem_2d_plots["p[e+] v p[e-]"] 
        = new TH2F("pep_v_pem", "pep_v_pem", 50, 0, 2.0, 50, 0, 2.0);
    
    track_emem_2d_plots["p[e-] v p[e-]"] 
        = new TH2F("pem_v_pem", "pem_v_pem", 50, 0, 2.0, 50, 0, 2.0);
    track_emem_2d_plots["theta[e-] v theta[e-]"] 
        = new TH2F("theta_em_v_theta_em", "theta_em_v_theta_em", 100, -0.05, 0.05, 100, -0.05, 0.05); 

}

std::string TrackAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;   
}

double TrackAnalysis::getMagnitude(std::vector<double> v) { 
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}
