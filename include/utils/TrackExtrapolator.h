
#ifndef __TRACK_EXTRAPOLATOR_H__
#define __TRACK_EXTRAPOLATOR_H__

#include <cmath>

#include <SvtTrack.h>

#include <AnalysisUtils.h>

namespace TrackExtrapolator { 

    /**
     *
     */
    double getR(SvtTrack* track);

    /**
     *
     */
    double getX0(SvtTrack* track);

    /**
     *
     */
    double getY0(SvtTrack* track);

    /**
     *
     */
    double getXc(SvtTrack* track);

    /**
     *
     */
    double getYc(SvtTrack* track);

    /**
     *
     */
    double getPathLength(SvtTrack* track, double x1, double y1, double x2, double y2);

    /**
     *
     */
    double getPathToXPlane(SvtTrack* track, double x);

    /**
     *
     */
    std::vector<double> getPointOnHelix(SvtTrack* track, double path_length);

    /**
     *
     */
    std::vector<double> extrapolateHelixToXPlane(SvtTrack* track, double x);    

}

#endif // __TRACK_EXTRAPOLATOR_H__
