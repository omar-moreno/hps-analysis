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
      canvas(NULL) { 
}

TrackAnalysis::~TrackAnalysis() { 

    delete canvas;

    std::unordered_map<std::string, TH1F*>::iterator track_parameter_plots_it = track_parameter_plots.begin();
    for (track_parameter_plots_it; 
            track_parameter_plots_it != track_parameter_plots.end(); 
            ++track_parameter_plots_it) { 
        delete track_parameter_plots_it->second;
    }
    track_parameter_plots.clear();
}

void TrackAnalysis::initialize() { 
}

void TrackAnalysis::processEvent(HpsEvent* event) {

    // Loop over all of the tracks in the event
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 

        // Get a track from the event
        track = event->getTrack(track_n); 

        track_parameter_plots["doca"]->Fill(track->getD0());
        track_parameter_plots["z0"]->Fill(track->getZ0());
        track_parameter_plots["phi0"]->Fill(track->getPhi());
        track_parameter_plots["curvature"]->Fill(track->getOmega());
        track_parameter_plots["tan_lambda"]->Fill(track->getTanLambda());
        track_parameter_plots["chi2"]->Fill(track->getChi2());
    }
}

void TrackAnalysis::finalize() { 
}

void TrackAnalysis::bookHistograms() { 

    canvas = new TCanvas("canvas", "canvas", 500, 500);
    
    track_parameter_plots["doca"] = new TH1F("doca", "doca", 80, -80, 80);
    track_parameter_plots["z0"] = new TH1F("z0", "z0", 30, -30, 30);
    track_parameter_plots["phi0"] = new TH1F("phi0", "phi0", 50, -0.5, 0.5);
    track_parameter_plots["curvature"] = new TH1F("curvature", "curvature", 200, -1, 1);
    track_parameter_plots["tan_lambda"] = new TH1F("tan_lambda", "tan_lambda", 100, -1, 1);
    track_parameter_plots["chi2"] = new TH1F("chi2", "chi2", 100, 0, 100);
}

std::string TrackAnalysis::toString() { 
}
