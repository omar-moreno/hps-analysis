
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

    if (z >= 997.2) { 
       position = extrapolateHelixToXPlane(track, 997.2);
       dz = z - 997.2;
    } else if (z <= 0) { 
       position = extrapolateHelixToXPlane(track, 0);
       dz = z - position[0];  
    } else { 
        position = extrapolateHelixToXPlane(track, z); 
        return position;
    }

    double phi = getPhi(track, position);

    double r = dz/getSinTheta(track)*cos(phi);
    double dx = r*getSinTheta(track)*sin(phi);
    double dy = r*getCosTheta(track);

    position[1] = position[1] + dx; 
    position[2] = position[2] + dy;

    return position;
}

