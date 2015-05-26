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
      class_name("TrackAnalysis") {  
}

TrackAnalysis::~TrackAnalysis() { 

    delete canvas;
    delete output_file;

    track_parameter_plots.clear();
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

        track_plots["doca"]->Fill(track->getD0());
        track_plots["z0"]->Fill(track->getZ0());
        track_plots["phi0"]->Fill(track->getPhi());
        track_plots["curvature"]->Fill(track->getOmega());
        track_plots["tan_lambda"]->Fill(track->getTanLambda());
        
        
        std::vector<double> p = track->getMomentum();
        double p_mag = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        double pt = sqrt(p[0]*p[0] + p[1]*p[1]);
        track_plots["p"]->Fill(p_mag);
        track_plots["pt"]->Fill(pt);
        track_plots["px"]->Fill(p[0]);
        track_plots["py"]->Fill(p[1]);
        track_plots["pz"]->Fill(p[2]);
        track_plots_2d["pz v px"]->Fill(p[2], p[0]);
        track_plots_2d["pz v py"]->Fill(p[2], p[1]);
        track_plots_2d["p v pt"]->Fill(p_mag, pt);

        if (track->getCharge() < 0) { 
            electron_track_plots["p"]->Fill(p_mag);
            electron_track_plots["pt"]->Fill(pt);
            electron_track_plots["px"]->Fill(p[0]);
            electron_track_plots["py"]->Fill(p[1]);
            electron_track_plots["pz"]->Fill(p[2]);
            electron_track_plots["chi2"]->Fill(track->getChi2());
            electron_track_plots["doca"]->Fill(track->getD0());
            electron_track_plots["z0"]->Fill(track->getZ0());
            electron_track_plots["phi0"]->Fill(track->getPhi());
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
            positron_track_plots["phi0"]->Fill(track->getPhi());
            positron_track_plots["curvature"]->Fill(track->getOmega());
            positron_track_plots["tan_lambda"]->Fill(track->getTanLambda());
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
    std::unordered_map<std::string, TH1F*>::iterator track_parameter_plots_it = track_parameter_plots.begin();
    for (track_parameter_plots_it; 
            track_parameter_plots_it != track_parameter_plots.end(); 
            ++track_parameter_plots_it) { 
        track_parameter_plots_it->second->Draw();
        track_parameter_plots_it->second->Write();
        canvas->Print("track_analysis_results.pdf("); 
    }

    std::unordered_map<std::string, TH1F*>::iterator track_plots_it = track_plots.begin();
    for (track_plots_it; 
            track_plots_it != track_plots.end(); 
            ++track_plots_it) { 
        track_plots_it->second->Draw();
        track_plots_it->second->Write();
        
        if (electron_track_plots[track_plots_it->first] != NULL) { 
            electron_track_plots[track_plots_it->first]->Draw("same");
            electron_track_plots[track_plots_it->first]->Write();
        }
        
        if (positron_track_plots[track_plots_it->first] != NULL) { 
            positron_track_plots[track_plots_it->first]->Draw("same");
            positron_track_plots[track_plots_it->first]->Write();
        }

        canvas->Print("track_analysis_results.pdf("); 
    }
    
    electron_track_plots["ep"]->Draw();
    canvas->Print("track_analysis_results.pdf("); 
    electron_track_plots["ep"]->Write();


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
    track_plots["phi0"] = new TH1F("phi0", "phi0", 40, -0.4, 0.4);
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
    electron_track_plots["phi0"] = new TH1F("phi0_electron", "phi0_electron", 40, -0.4, 0.4);
    electron_track_plots["curvature"] = new TH1F("curvature_electron", "curvature_electron", 25, -0.0025, 0.0025);
    electron_track_plots["tan_lambda"] = new TH1F("tan_lambda_electron", "tan_lambda_electron", 100, -0.1, 0.1);

    electron_track_plots["p"] = new TH1F("p_electron", "p_electron", 50, 0, 2.0);
    electron_track_plots["pt"] = new TH1F("pt_electron", "pt_electron", 50, -0.1, 0.2);
    electron_track_plots["px"] = new TH1F("px_electron", "px_electron", 50, -0.1, 0.2); 
    electron_track_plots["py"] = new TH1F("py_electron", "py_electron", 50, -0.2, 0.2); 
    electron_track_plots["pz"] = new TH1F("pz_electron", "pz_electron", 50, 0, 3.0); 
    electron_track_plots["chi2"] = new TH1F("chi2_electron", "chi2_electron", 40, 0, 40);
    electron_track_plots["ep"] = new TH1F("ep", "ep", 60, 0, 2);

    positron_track_plots["doca"] = new TH1F("doca_positron", "doca_positron", 40, -20, 20);
    positron_track_plots["z0"] = new TH1F("z0_positron", "z0_positron", 100, -5, 5);
    positron_track_plots["phi0"] = new TH1F("phi0_positron", "phi0_positron", 40, -0.4, 0.4);
    positron_track_plots["curvature"] = new TH1F("curvature_positron", "curvature_positron", 25, -0.0025, 0.0025);
    positron_track_plots["tan_lambda"] = new TH1F("tan_lambda_positron", "tan_lambda_positron", 100, -0.1, 0.1);

    positron_track_plots["p"] = new TH1F("p_positron", "p_positron", 50, 0, 2.0);
    positron_track_plots["pt"] = new TH1F("pt_positron", "pt_positron", 50, -0.1, 0.2);
    positron_track_plots["px"] = new TH1F("px_positron", "px_positron", 50, -0.1, 0.2); 
    positron_track_plots["py"] = new TH1F("py_positron", "py_positron", 50, -0.15, 0.15); 
    positron_track_plots["pz"] = new TH1F("pz_positron", "pz_positron", 50, 0, 3.0); 
    positron_track_plots["chi2"] = new TH1F("chi2_positron", "chi2_positron", 40, 0, 40);

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
