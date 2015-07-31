/**
 * @file AnalysisUtils.h
 * @brief A set of utilities commonly used when doing analysis
 * @author Omar Moreno <omoreno1@ucsc.edu> \n
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date May 20, 2015
 */

#ifndef __ANALYSIS_UTILS_H__
#define __ANALYSIS_UTILS_H__

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <vector>
#include <cmath>

//---------------//
//--- HPS DST ---//
//---------------//
#include <HpsEvent.h>
#include <EcalCluster.h>
#include <SvtTrack.h>

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <TrackExtrapolator.h>

namespace AnalysisUtils { 

    /**
     * Get the magnitude of a 3D vector.  If the vector doesn't
     * contain exactly three elements, an exception is thrown.
     *
     * @param v : 3D vector
     * @return Magnitude of the vector
     */
    double getMagnitude(std::vector<double> v);

    /**
     * Calculate the invariant mass of a pair of tracks.
     *
     * @param track_0 SVT track composing a pair
     * @param track_1 SVT track composing a pair
     * @return invariant mass of the track pair
     */
    double getInvariantMass(SvtTrack* track_0, SvtTrack* track_1);

    /**
     * Get a cluster pair from an event.
     *
     * @param event An HpsEvent object
     * @return a vector containing the pair of clusters
     */    
    std::vector<EcalCluster*> getClusterPair(HpsEvent* event);

    /**
     *
     */
    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

}

#endif // __ANALYSIS_UTILS_H__
