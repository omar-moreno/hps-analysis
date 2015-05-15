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

    /*std::unordered_map<std::string, TH1F*>::iterator track_parameter_plots_it = track_parameter_plots.begin();
    for (track_parameter_plots_it; 
            track_parameter_plots_it != track_parameter_plots.end(); 
            ++track_parameter_plots_it) { 
        delete track_parameter_plots_it->second;
    }*/
    track_parameter_plots.clear();
}

void TrackAnalysis::initialize() { 

    this->bookHistograms();
}

void TrackAnalysis::processEvent(HpsEvent* event) {

    track_plots["Number of tracks"]->Fill(event->getNumberOfTracks());

    // Loop over all of the tracks in the event
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 

        // Get a track from the event
        track = event->getTrack(track_n); 

        track_parameter_plots["doca"]->Fill(track->getD0());
        track_parameter_plots["z0"]->Fill(track->getZ0());
        track_parameter_plots["phi0"]->Fill(track->getPhi());
        track_parameter_plots["curvature"]->Fill(track->getOmega());
        track_parameter_plots["tan_lambda"]->Fill(track->getTanLambda());
        
        track_plots["chi2"]->Fill(track->getChi2());
    
        track_plots["Hits per track"]->Fill(track->getSvtHits()->GetEntriesFast()); 
    }
}

void TrackAnalysis::finalize() {

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
        canvas->Print("track_analysis_results.pdf("); 
    }
    
    canvas->Print("track_analysis_results.pdf]");
    output_file->Close(); 
}

void TrackAnalysis::bookHistograms() { 

    canvas = new TCanvas("canvas", "canvas", 500, 500);
    
    track_parameter_plots["doca"] = new TH1F("doca", "doca", 80, -80, 80);
    track_parameter_plots["z0"] = new TH1F("z0", "z0", 30, -30, 30);
    track_parameter_plots["phi0"] = new TH1F("phi0", "phi0", 50, -0.5, 0.5);
    track_parameter_plots["curvature"] = new TH1F("curvature", "curvature", 200, -1, 1);
    track_parameter_plots["tan_lambda"] = new TH1F("tan_lambda", "tan_lambda", 100, -1, 1);

    track_plots["Number of tracks"] = new TH1F("number_of_tracks", "number_of_tracks", 10, 0, 10);
    track_plots["Hits per track"] = new TH1F("hits_per_track", "hits_per_track", 6, 1, 7);
    track_plots["Track charge"] = new TH1F("track_charge", "track_charge", 3, -1, 2);
    track_plots["chi2"] = new TH1F("chi2", "chi2", 100, 0, 100);
    track_plots["px"] = new TH1F("px", "px", 50, -0.1, 0.2); 
    track_plots["py"] = new TH1F("py", "py", 50, -0.2, 0.2); 
    track_plots["pz"] = new TH1F("pz", "pz", 50, 0, 3.0); 
}

std::string TrackAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;   
}
