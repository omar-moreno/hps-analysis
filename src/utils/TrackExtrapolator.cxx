
#include <TrackExtrapolator.h>

double TrackExtrapolator::getR(SvtTrack* track) {
    return 1/track->getOmega();
}

double TrackExtrapolator::getX0(SvtTrack* track) {
    return -(track->getD0()*sin(track->getPhi0()));
}

double TrackExtrapolator::getY0(SvtTrack* track) {
    return track->getD0()*cos(track->getPhi0());
}

double TrackExtrapolator::getXc(SvtTrack* track) { 
    return (getR(track) - track->getD0())*sin(track->getPhi0());
}

double TrackExtrapolator::getYc(SvtTrack* track) { 
    return -(getR(track) - track->getD0())*cos(track->getPhi0());
}

double TrackExtrapolator::getPathLength(SvtTrack* track, double x1, double y1, double x2, double y2) { 
    double phi1 = atan2(y1 - getYc(track), x1 - getXc(track));
    double phi2 = atan2(y2 - getYc(track), x2 - getXc(track));
    double dphi = phi2 - phi1;

    if (dphi > 3.14159) dphi -= 2*3.14159;
    else if (dphi < -3.14159) dphi += 2*3.14159;

    return -getR(track)*dphi;
}

double TrackExtrapolator::getPathToXPlane(SvtTrack* track, double x) {

    double r = getR(track);

    double y = getYc(track) + AnalysisUtils::sgn<double>(r)*sqrt(r*r - pow(x - getXc(track), 2));

    return getPathLength(track, getX0(track), getY0(track), x, y); 
}

double TrackExtrapolator::getPhi(SvtTrack* track, std::vector<double> position) { 

    double x = sin(track->getPhi0()) - track->getOmega()*(position[0] - getX0(track));
    double y = cos(track->getPhi0()) + track->getOmega()*(position[1] - getY0(track));

    return atan2(x, y);
}

double TrackExtrapolator::getSinTheta(SvtTrack* track) { 
    return 1/sqrt(1 + pow(track->getTanLambda(), 2));
}

double TrackExtrapolator::getCosTheta(SvtTrack* track) { 
    return track->getTanLambda()/sqrt(1 + pow(track->getTanLambda(), 2));
}

std::vector<double> TrackExtrapolator::getPointOnHelix(SvtTrack* track, double path_length) {

    double phi = track->getPhi0() - (path_length/getR(track));
    double x = getXc(track) - getR(track)*sin(phi); 
    double y = getYc(track) + getR(track)*cos(phi);
    double z = track->getZ0() + path_length*track->getTanLambda();

    std::vector<double> position(3, 0);
    position[0] = x; 
    position[1] = y;
    position[2] = z;

    return position;  
}

std::vector<double> TrackExtrapolator::extrapolateHelixToXPlane(SvtTrack* track, double x) { 

    double path_length = getPathToXPlane(track, x);

    return getPointOnHelix(track, path_length);
}

std::vector<double> TrackExtrapolator::extrapolateTrack(SvtTrack* track, double z) { 

    std::vector<double> position(3,0); 
    double dz = 0;

    if (z >= TrackExtrapolator::DIPOLE_EDGE) {

        // If the point of extrapolation is outside of the dipole edge, then 
        // extrapolate the helix to the edge and then use a straight line 
        // extrapolation beyond that

        // Extrapolate the helix to the edge of the dipole 
        position = extrapolateHelixToXPlane(track, TrackExtrapolator::DIPOLE_EDGE);

        // Get the difference between the dipole edge and the extrapolation
        // point. The track will be extrapolated assuming no field for this
        // distance i.e. straight line extrapolation
        dz = z - TrackExtrapolator::DIPOLE_EDGE;
    } else if (z <= 0) {

        // If the extrapolation point is upstream of the target, do something
        // similar as above

        position = extrapolateHelixToXPlane(track, 0);
        dz = z - position[0];  
    
    } else { 
    
        // If the extrapolation point is inside of the field region, 
        // analytically extrapolate the helix and return the position
        position = extrapolateHelixToXPlane(track, z); 
   
        // FIXME: This position should be in the JLab frame     
        return position;
    }

    // Calculate the value of Phi at the track position
    double phi = getPhi(track, position);
    
    // Calcualte the distance to the extrapolation point
    double r = dz/getSinTheta(track)*cos(phi);

    // Get the delta x and y values at the point of extrapolation 
    double dx = r*getSinTheta(track)*sin(phi);
    double dy = r*getCosTheta(track);

    // Calculate the position of the track at the extrapolation point
    std::vector<double> extrapolated_position(3, 0);
    extrapolated_position[0] = position[1] + dx;
    extrapolated_position[1] = position[2] + dy;
    extrapolated_position[2] = z; 

    return extrapolated_position;
}
