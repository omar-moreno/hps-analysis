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

    // Loop over all of the tracks in the event
    for (int first_track_n = 0; first_track_n < event->getNumberOfTracks(); ++first_track_n) { 

        // Get a track from the event
        first_track = event->getTrack(first_track_n);

        // Get the stereo hits associated with the track
        TRefArray* first_track_hits = first_track->getSvtHits();

        for (int second_track_n = 0; second_track_n < event->getNumberOfTracks(); ++second_track_n) { 
           
            second_track = event->getTrack(second_track_n);

            if (first_track == second_track) { 
                //std::cout << "Cannot compare hits from the same track." << std::endl;
                continue;
            }

            TRefArray* second_track_hits = second_track->getSvtHits(); 
            
            for (int hit_n = 0; hit_n < first_track_hits->GetSize(); ++hit_n) { 
                
                SvtHit* first_track_hit = (SvtHit*) first_track_hits->At(hit_n);

                for (int second_track_hit_n = 0; second_track_hit_n < second_track_hits->GetSize(); ++second_track_hit_n) {
                    
                    SvtHit* second_track_hit = (SvtHit*) second_track_hits->At(second_track_hit_n);

                    if (first_track_hit->getPosition()[2] == second_track_hit->getPosition()[2]) {
                        std::cout << "Shared Hit Found" << std::endl;
                    }
                }
            }
        } 
    }
}

void SharedHitAnalysis::finalize() { 

}

void SharedHitAnalysis::bookHistograms() { 

}

std::string SharedHitAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;   
}
