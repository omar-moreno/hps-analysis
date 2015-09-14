/**
 *
 * @file TrackType.cxx
 * @brief A set of utilities to determine what track finding algorithm was
 *        used to find a track and whether it was refit using GBL.
 * @author <a href="mailto:omoreno1@ucsc.edu">Omar Moreno</a>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date September 13, 2015
 *
 */

#include <TrackType.h>

bool TrackType::foundByStrategy(SvtTrack* track, StrategyType strategyType) { 
    return (track->getType() & (1 << ( ((int) strategyType) - 1))) != 0; 
}

bool TrackType::isGbl(SvtTrack* track) { 
    return (track->getType() & (1 << 5)) != 0; 
}
