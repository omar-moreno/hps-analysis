/**
 *
 * @file TrackUtils.cxx
 * @brief A set of {@link SvtTrack} utilities. 
 * @author <a href="mailto:omoreno1@ucsc.edu">Omar Moreno</a>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date September 15, 2015
 *
 */

#include <TrackUtils.h>

std::vector<SvtTrack*> TrackUtils::getGoodTracksList(HpsEvent* event) { 

    //std::cout << "Event: " << event->getEventNumber() << std::endl;
    //std::cout << "Get number of tracks: " << event->getNumberOfTracks() << std::endl;

    // Map the hits on a track to the corresponding track
    std::map<std::vector<SvtHit*>, std::vector<SvtTrack*>> hit_to_tracks_map;

    std::vector<SvtTrack*> tracks;

    // Loop over all of the tracks and remove GBL tracks and duplicates
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 

        // Get a track from the event 
        SvtTrack* track = event->getTrack(track_n);

        // Only add seed tracks to the list of 'good' tracks
        if (TrackType::isGbl(track)) { 
            //std::cout << "GBL track found! Skipping it" << std::endl;
            continue; 
        }
    
        // Get the hits on a track
        TRefArray* track_hits = track->getSvtHits(); 

        // Add the hits to a std::vector.  This is done to take advantage of 
        // a std::vector's == operator.
        std::vector<SvtHit*> hit_vector(6, nullptr);  
        for (int hit_n = 0; hit_n < track_hits->GetSize(); ++hit_n) { 
            SvtHit* hit = (SvtHit*) track_hits->At(hit_n);
            hit_vector.at(hit->getLayer()-1) = hit;   
        }

        // Add the list of hits to the map if it doesn't exist already
        auto it = hit_to_tracks_map.find(hit_vector); 
        if (it == hit_to_tracks_map.end()) { 
            //std::cout << "Adding hits to map: " << std::endl;
            std::vector<SvtTrack*> track_list; 
            hit_to_tracks_map.emplace(hit_vector, track_list);  
        }

        // Add the track to the map
        hit_to_tracks_map[hit_vector].push_back(track);
        //std::cout << "Track size: " << hit_to_tracks_map[hit_vector].size() << std::endl; 
    }

    // Keep the best track of the duplicates. For now, best is defined as the 
    // track with the best chi2. 
    //std::cout << "Comparing equivalent tracks: " << std::endl;
    for (auto& map_item : hit_to_tracks_map) { 
        int index = 0; 
        double chi2 = 1000;
        SvtTrack* best_track = nullptr;  
        for (auto& duplicate_track : map_item.second) { 
            if (duplicate_track->getChi2() < chi2) { 
                chi2 = duplicate_track->getChi2(); 
                best_track = duplicate_track; 
            }
            //std::cout << "Track: " << index << " Chi2 " << duplicate_track->getChi2() << std::endl;
        }
        //std::cout << "Chi2 of best track: " << best_track->getChi2() << std::endl;
        tracks.push_back(best_track); 
    } 

    return tracks; 
}
