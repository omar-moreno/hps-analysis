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

    // Loop over all of the tracks and remove GBL tracks and duplicates
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 

        // Get a track from the event 
        SvtTrack* track = event->getTrack(track_n);

        // Only add seed tracks to the list of 'good' tracks
        if (TrackType::isGbl(track)) continue; 
    
        // Get the hits associated with a track
        TRefArray* track_hits = track->getSvtHits(); 

        // Add the hits to a std::vector.  This is done to take advantage of 
        // a std::vector's == operator. A TRefArray doesn't have a define
        // way to compare references to objects
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

    // Container to hold the 'good' tracks
    std::vector<SvtTrack*> tracks;

    // Loop through the map of hits to tracks and keep the best track of those
    // which are duplicates i.e. they share the same hits. For now, best is
    // defined as the track with the best chi2. 
    //std::cout << "Comparing equivalent tracks: " << std::endl;
    for (auto& map_item : hit_to_tracks_map) { 
       
        if (map_item.second.size() == 1) { 
            tracks.push_back(map_item.second[0]); 
            continue; 
        }

        SvtTrack* best_track = map_item.second[0];  
        for (auto& track : map_item.second) { 
            if (track->getChi2() < best_track->getChi2()) { 
                best_track = track; 
            }
        }
        tracks.push_back(best_track); 
    } 

    // Loop through all of the remaining 'good' tracks and remove 
    // 'partial tracks' i.e. tracks whose hit content is a subset of another
    // track.
    for (auto& track : tracks) { 
        
    }

    return tracks; 
}
