/**
 * @file SharedHitAnalysis.h
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date June 15, 2015
 */

#include <SharedHitAnalysis.h>

SharedHitAnalysis::SharedHitAnalysis()
    : track_plotter(new Plotter()),
      electron_plotter(new Plotter()), 
      positron_plotter(new Plotter()), 
      class_name("SharedHitAnalysis") { 
    
}

SharedHitAnalysis::~SharedHitAnalysis() { 
    delete track_plotter;
}

void SharedHitAnalysis::initialize() { 
    this->bookHistograms();
}

void SharedHitAnalysis::processEvent(HpsEvent* event) {

    // Copy the tracks to a list
    std::list<SvtTrack*> tracks;  
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) {
         tracks.push_back(event->getTrack(track_n));
    }
    //std::cout << "Total number of tracks: " << tracks.size() << std::endl;

    // Loop through all the tracks and find all tracks that share a hit
    std::list<SvtTrack*>::iterator track_it = tracks.begin();
    int first_track_index = 0;
    std::list<SvtHit*> shared_hits;
    while (track_it != tracks.end()) { 
        
        //std::cout << "First track index: " << ++first_track_index << std::endl;

        //std::list<SvtTrack*> shared_hit_tracks; 
        
        // Get the stereo hits associated with the track
        TRefArray* track_hits = (*track_it)->getSvtHits();

        int second_track_index = 0; 
        std::list<SvtTrack*>::iterator track_cand_it = tracks.begin();
        for (track_cand_it; track_cand_it != tracks.end(); ++track_cand_it) { 
            
            //std::cout << "Second track index: " << ++second_track_index << std::endl;
            if (track_it == track_cand_it) { 
                //std::cout << "Cannot compare hits from the same track." << std::endl;
                continue;
            }
            
            TRefArray* track_cand_hits = (*track_cand_it)->getSvtHits();
            // Loop through the hits of the first track and compare them to the hits of 
            // the second track
            int shared_hit_count = 0;
            for (int hit_n = 0; hit_n < track_hits->GetSize(); ++hit_n) {
                
                SvtHit* track_hit = (SvtHit*) track_hits->At(hit_n);
                
                for (int cand_hit_n = 0; cand_hit_n < track_cand_hits->GetSize(); ++cand_hit_n) { 
                    
                    SvtHit* track_cand_hit = (SvtHit*) track_cand_hits->At(cand_hit_n);     
                       
                    if (track_hit->getPosition()[2] == track_cand_hit->getPosition()[2]) {
                        
                        //std::cout << "Shared hit found on layer " << track_cand_hit->getLayer() << std::endl;
                        track_plotter->get1DHistogram("shared hit layer")->Fill(track_hit->getLayer());
                        shared_hit_count++;
                        if (std::find(shared_hits.begin(), shared_hits.end(), track_hit) == shared_hits.end())  
                            shared_hits.push_back(track_hit);
                    }
                }
            }
            if (shared_hit_count > 0) { 
                track_plotter->get1DHistogram("shared hits between tracks")->Fill(shared_hit_count);
            }
        }
        tracks.erase(track_it++);
    }

    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) {
        
        int shared_hit_count = 0;
        SvtTrack* track = event->getTrack(track_n);

        // Get the stereo hits associated with the track
        TRefArray* track_hits = track->getSvtHits();
        for (int hit_n = 0; hit_n < track_hits->GetSize(); ++hit_n) {
                
            SvtHit* track_hit = (SvtHit*) track_hits->At(hit_n);
            
            if (std::find(shared_hits.begin(), shared_hits.end(), track_hit) != shared_hits.end()) {
                shared_hit_count++; 
            }  
        }
   
        track_plotter->get1DHistogram("doca - shared " + std::to_string(shared_hit_count))->Fill(track->getD0());
        track_plotter->get1DHistogram("z0 - shared " + std::to_string(shared_hit_count))->Fill(track->getZ0());
        track_plotter->get1DHistogram("sin(phi0) - shared " + std::to_string(shared_hit_count))->Fill(sin(track->getPhi0()));
        track_plotter->get1DHistogram("curvature - shared " + std::to_string(shared_hit_count))->Fill(track->getOmega());
        track_plotter->get1DHistogram("tan_lambda - shared " + std::to_string(shared_hit_count))->Fill(track->getTanLambda());
           
        track_plotter->get1DHistogram("shared hits per track")->Fill(shared_hit_count);

        if (track->getCharge() < 0) { 
            
            electron_plotter->get1DHistogram("doca - shared " + std::to_string(shared_hit_count))->Fill(track->getD0());
            electron_plotter->get1DHistogram("z0 - shared " + std::to_string(shared_hit_count))->Fill(track->getZ0());
            electron_plotter->get1DHistogram("sin(phi0) - shared " + std::to_string(shared_hit_count))->Fill(sin(track->getPhi0()));
            electron_plotter->get1DHistogram("curvature - shared " + std::to_string(shared_hit_count))->Fill(track->getOmega());
            electron_plotter->get1DHistogram("tan_lambda - shared " + std::to_string(shared_hit_count))->Fill(track->getTanLambda());
        } else { 
            positron_plotter->get1DHistogram("doca - shared " + std::to_string(shared_hit_count))->Fill(track->getD0());
            positron_plotter->get1DHistogram("z0 - shared " + std::to_string(shared_hit_count))->Fill(track->getZ0());
            positron_plotter->get1DHistogram("sin(phi0) - shared " + std::to_string(shared_hit_count))->Fill(sin(track->getPhi0()));
            positron_plotter->get1DHistogram("curvature - shared " + std::to_string(shared_hit_count))->Fill(track->getOmega());
            positron_plotter->get1DHistogram("tan_lambda - shared " + std::to_string(shared_hit_count))->Fill(track->getTanLambda());
        }
    }
}

void SharedHitAnalysis::finalize() { 
    track_plotter->saveToPdf("shared_hit_analysis.pdf");
    track_plotter->saveToRootFile("shared_hit_analysis.root");

    electron_plotter->saveToPdf("shared_hit_analysis_electron.pdf");
    electron_plotter->saveToRootFile("shared_hit_analysis_electron.root");

    positron_plotter->saveToPdf("shared_hit_analysis_positron.pdf");
    positron_plotter->saveToRootFile("shared_hit_analysis_positron.root");
}

void SharedHitAnalysis::bookHistograms() { 

    track_plotter->build1DHistogram("shared hits between tracks", 6, 0, 6);
    track_plotter->build1DHistogram("shared hits per track", 6, 0, 6);
    track_plotter->build1DHistogram("shared hit layer", 6, 1, 7);
    //track_plotter->build1DHistogram("first track #chi^{2} - second track #chi^{2}", 100, -25, 25);
    //track_plotter->build2DHistogram("first track v second track #chi^{2}", 50, 0, 25, 50, 0, 25);

    for (int shared_hit_n = 0; shared_hit_n < 7; ++shared_hit_n) { 
       
        track_plotter->build1DHistogram("doca - shared " + std::to_string(shared_hit_n), 80, -10, 10);
        track_plotter->build1DHistogram("z0 - shared " + std::to_string(shared_hit_n), 80, -2, 2);
        track_plotter->build1DHistogram("sin(phi0) - shared " + std::to_string(shared_hit_n), 40, -0.2, 0.2);
        track_plotter->build1DHistogram("curvature - shared " + std::to_string(shared_hit_n), 50, -0.001, 0.001);
        track_plotter->build1DHistogram("tan_lambda - shared " + std::to_string(shared_hit_n), 100, -0.1, 0.1);
        
        electron_plotter->build1DHistogram("doca - shared " + std::to_string(shared_hit_n), 80, -10, 10);
        electron_plotter->build1DHistogram("z0 - shared " + std::to_string(shared_hit_n), 80, -2, 2);
        electron_plotter->build1DHistogram("sin(phi0) - shared " + std::to_string(shared_hit_n), 40, -0.2, 0.2);
        electron_plotter->build1DHistogram("curvature - shared " + std::to_string(shared_hit_n), 50, -0.001, 0.001);
        electron_plotter->build1DHistogram("tan_lambda - shared " + std::to_string(shared_hit_n), 100, -0.1, 0.1);
        
        positron_plotter->build1DHistogram("doca - shared " + std::to_string(shared_hit_n), 80, -10, 10);
        positron_plotter->build1DHistogram("z0 - shared " + std::to_string(shared_hit_n), 80, -2, 2);
        positron_plotter->build1DHistogram("sin(phi0) - shared " + std::to_string(shared_hit_n), 40, -0.2, 0.2);
        positron_plotter->build1DHistogram("curvature - shared " + std::to_string(shared_hit_n), 50, -0.001, 0.001);
        positron_plotter->build1DHistogram("tan_lambda - shared " + std::to_string(shared_hit_n), 100, -0.1, 0.1);

    } 
}

std::string SharedHitAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;   
}
