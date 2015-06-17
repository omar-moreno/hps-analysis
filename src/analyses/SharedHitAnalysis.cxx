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
    : first_track(NULL),
      second_track(NULL),
      track_plotter(new Plotter()), 
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
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) {
    
    }


    // Loop over all of the tracks in the event
    for (int first_track_n = 0; first_track_n < event->getNumberOfTracks(); ++first_track_n) { 

        // Get a track from the event
        first_track = event->getTrack(first_track_n);

        //std::cout << "First track chi^2: " << first_track->getChi2() << std::endl;

        // Get the stereo hits associated with the track
        TRefArray* first_track_hits = first_track->getSvtHits();

        for (int second_track_n = 0; second_track_n < event->getNumberOfTracks(); ++second_track_n) { 
           
            second_track = event->getTrack(second_track_n);

            if (first_track == second_track) { 
                //std::cout << "Cannot compare hits from the same track." << std::endl;
                continue;
            }

            if (std::find(checked_tracks.begin(), checked_tracks.end(), second_track) != checked_tracks.end()) continue;

            //std::cout << "Second track chi^2: " << second_track->getChi2() << std::endl;
            TRefArray* second_track_hits = second_track->getSvtHits(); 
           
            double shared_hit = 0; 
            for (int hit_n = 0; hit_n < first_track_hits->GetSize(); ++hit_n) { 
                
                SvtHit* first_track_hit = (SvtHit*) first_track_hits->At(hit_n);

                for (int second_track_hit_n = 0; second_track_hit_n < second_track_hits->GetSize(); ++second_track_hit_n) {
                    
                    SvtHit* second_track_hit = (SvtHit*) second_track_hits->At(second_track_hit_n);

                    if (first_track_hit->getPosition()[2] == second_track_hit->getPosition()[2]) {
                        track_plotter->get1DHistogram("shared hit layer")->Fill(second_track_hit->getLayer());
                        //std::cout << "Shared hit found on layer " << second_track_hit->getLayer() << std::endl;
                        shared_hit++;
                    }
                }
            }
            if (shared_hit > 0) { 
                track_plotter->get1DHistogram("shared hits per track")->Fill(shared_hit); 
                track_plotter->get2DHistogram("first track v second track #chi^{2}")->Fill(first_track->getChi2(), second_track->getChi2());
                track_plotter->get1DHistogram("first track #chi^{2} - second track #chi^{2}")->Fill(first_track->getChi2() - second_track->getChi2());
            } 
        } 
        checked_tracks.push_back(first_track);
    }
}

void SharedHitAnalysis::finalize() { 
    track_plotter->saveToPdf("shared_hit_analysis.pdf");
}

void SharedHitAnalysis::bookHistograms() { 

    track_plotter->build1DHistogram("shared hits per track", 6, 0, 6);
    track_plotter->build1DHistogram("shared hit layer", 6, 1, 7);
    track_plotter->build1DHistogram("first track #chi^{2} - second track #chi^{2}", 100, -25, 25);
    track_plotter->build2DHistogram("first track v second track #chi^{2}", 50, 0, 25, 50, 0, 25);
}

std::string SharedHitAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;   
}
