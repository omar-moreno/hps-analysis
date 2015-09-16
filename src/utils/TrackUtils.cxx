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
    
    std::cout << "Event: " << event->getEventNumber() << std::endl;
    std::cout << "Get number of tracks: " << event->getNumberOfTracks() << std::endl;

    std::vector<SvtTrack*> tracks; 
    std::map<std::vector<SvtHit*>, std::vector<SvtTrack*>> equivalent_tracks;

    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
        
        SvtTrack* track = event->getTrack(track_n);
        
        if (TrackType::isGbl(track)) { 
            std::cout << "GBL track found! Skipping it" << std::endl;
            continue; 
        }
        
        TRefArray* track_hits = track->getSvtHits(); 
        std::vector<SvtHit*> hits(6, nullptr);  
        for (int hit_n = 0; hit_n < track_hits->GetSize(); ++hit_n) { 
            SvtHit* hit = (SvtHit*) track_hits->At(hit_n);
            hits.at(hit->getLayer()-1) = hit;   
        }

        auto it = equivalent_tracks.find(hits); 
        if (it == equivalent_tracks.end()) { 
            std::cout << "Adding hits to map: " << std::endl;
            std::vector<SvtTrack*> track_list; 
            equivalent_tracks.emplace(hits, track_list);  
        }

        equivalent_tracks[hits].push_back(track);
        std::cout << "Track size: " << equivalent_tracks[hits].size() << std::endl; 
    }

    std::cout << "Comparing equivalent tracks: " << std::endl;
    for (auto& map_item : equivalent_tracks) { 
        int index = 0; 
        for (auto& equivalent_track : map_item.second) { 
            std::cout << "Track: " << index << " Chi2 " << equivalent_track->getChi2() << std::endl;
        }
    } 

    return tracks; 
}
