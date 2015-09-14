/**
 *
 * @file TrackType.h
 * @brief A set of utilities to determine what track finding algorithm was
 *        used to find a track and whether it was refit using GBL.
 * @author <a href="mailto:omoreno1@ucsc.edu">Omar Moreno</a>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date September 13, 2015
 *
 */

#ifndef __TRACK_TYPE_H__
#define __TRACK_TYPE_H__

//-------------//
//   HPS DST   //
//-------------//
#include <SvtTrack.h>

//------------------//
//   HPS Analysis   //
//------------------//   
#include <StrategyType.h>

namespace TrackType { 
    
    /**
     * Checks the track type of an {@link SvtTrack} to determine if it was 
     * found by the given strategy.
     *
     * @param track {@link SvtTrack} whose track type will be checked.
     * @param strategyType {@link StrategyType} enum constant identifying the 
     *                     desired track finding algorithm.
     * @return Returns <code>true</code> if an {@link SvtTrack} was found by 
     *         the given strategy, <code>false</code> otherwise.
     */
    bool foundByStrategy(SvtTrack* track, StrategyType strategyType); 

    /**
     * Checks the track type of an {@link SvtTrack} to determine if it was 
     * refit using GBL.
     *
     * @param track {@link SvtTrack} whose track type will be checked.
     * @return Returns <code>true</code> if an {@link SvtTrack} was refit using
     *         GBL strategy, <code>false</code> otherwise.
     */
    bool isGbl(SvtTrack* track); 

}

#endif // __TRACK_TYPE_H__
