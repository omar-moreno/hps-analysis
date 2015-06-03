
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

