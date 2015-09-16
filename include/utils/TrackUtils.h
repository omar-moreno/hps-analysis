/**
 *
 * @file TrackUtils.h
 * @brief A set of {@link SvtTrack} utilities. 
 * @author <a href="mailto:omoreno1@ucsc.edu">Omar Moreno</a>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date September 15, 2015
 *
 */

#ifndef __TRACK_UTILS_H__
#define __TRACK_UTILS_H__

#include <vector>
#include <map>

#include <HpsEvent.h>
#include <SvtTrack.h>

#include <TrackType.h>

namespace TrackUtils { 

    /**
     * Get a 'good' list of tracks from the event.
     *
     * @param event The {@link HpsEvent} object from which the tracks will be 
     *              retrieved.
     * @return A vector containing references to 'good' {@link SvtTracks}.
     */
    std::vector<SvtTrack*> getGoodTracksList(HpsEvent* event); 

}

#endif // __TRACK_UTILS_H__
