
#ifndef __TRACK_EXTRAPOLATOR_H__
#define __TRACK_EXTRAPOLATOR_H__

#include <iostream>
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
        double getPhi(SvtTrack* track, std::vector<double> position);

        /**
         *
         */
        // TODO: Move this function to a track utility namespace
        double getSinTheta(SvtTrack* track); 

        /**
         *
         */
        // TODO: Move this function to a track utility namespace
        double getCosTheta(SvtTrack* track); 

        /**
         *
         */
        std::vector<double> getPointOnHelix(SvtTrack* track, double path_length);

        /**
         *
         */
        std::vector<double> extrapolateHelixToXPlane(SvtTrack* track, double x);   

        /**
         *
         */
        std::vector<double> extrapolateTrack(SvtTrack* track, double z);  

        /**
         * The end of the box dipole field during the engineering run.  Note 
         * that the dipole length was changed such that the integral of the 
         * box dipole field matches that of the field map integral.  The 
         * dipole edge is set to 45.72 cm (center of dipole) + 
         * 108/2 cm (dipole length).
         *
         */
        static const int DIPOLE_EDGE = 997.2;

}

#endif // __TRACK_EXTRAPOLATOR_H__
