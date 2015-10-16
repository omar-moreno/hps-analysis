'''
    A module used to extrapolate a HPS SVT track to any position within and 
    outside the analyzing magnet.

    Author : Omar Moreno <omoreno1@ucsc.edu>
    Institution : Santa Cruz Institute for Particle Physics
                  University of California, Santa Cruz
    Date : October 15, 2015

'''

import math

DIPOLE_EDGE = 997.2 # cm
DIPOLE_EDGE_LOW = 0.0 # cm

def get_R(track) : return 1/track.getOmega()

def get_xc(track) : return (get_R(track) - track.getD0())*math.sin(track.getPhi0())

def get_yc(track) : return -(get_R(track) - track.getD0())*math.cos(track.getPhi0())

def get_x0(track) : return -track.getD0()*math.sin(track.getPhi0())

def get_y0(track) : return track.getD0()*math.cos(track.getPhi0())

def get_path_length(track, x1, y1, x2, y2):
        
    phi1 = math.atan2(y1 - get_yc(track), x1 - get_xc(track))
    phi2 = math.atan2(y2 - get_yc(track), x2 - get_xc(track))
    dphi = phi2 - phi1

    if dphi > math.pi: dphi -= 2.*math.pi
    if dphi < -math.pi: dphi += 2.*math.pi

    return -get_R(track)*dphi

def get_path_to_x_plane(track, x):
        
    r = get_R(track)
    y = get_yc(track) + math.copysign(1, r)*math.sqrt(r*r - math.pow(x - get_xc(track), 2))

    return get_path_length(track, get_x0(track), get_y0(track), x, y)

def get_point_on_helix(track, path_length):
        
    r = get_R(track)
    phi = track.getPhi0() - (path_length/r)

    x = get_xc(track) - r*math.sin(phi)
    y = get_yc(track) + r*math.cos(phi)
    z = track.getZ0() + path_length*track.getTanLambda()

    return x, y, z

def extrapolate_helix_to_x_plane(track, x): 
        
    path_length = get_path_to_x_plane(track, x)
    return get_point_on_helix(track, path_length)


def get_phi(track, position):

	x = math.sin(track.getPhi0()) - track.getOmega()*(position[0] - get_x0(track))
	y = math.cos(track.getPhi0()) + track.getOmega()*(position[1] - get_y0(track))

	return math.atan2(x, y)

def get_sin_theta(track): return 1/math.sqrt(1 + math.pow(track.getTanLambda(),2))

def get_cos_theta(track) : return track.getTanLambda()/math.sqrt(1 + math.pow(track.getTanLambda(),2))

def extrapolate_track(track, z):

    if z >= DIPOLE_EDGE :

        '''
            If the point of extrapolation is outside of the dipole edge, then
            extrapolate the helix to the edge and then use a straight line 
            extrapolation beyond that.
        '''

        # Extrapolate the helix to the edge of the dipole
        position = extrapolate_helix_to_x_plane(track, DIPOLE_EDGE)

        '''
            Get the difference between the dipole edge and the extrapolation
            point.  The track will be extrapolated assuming no field for this
            distance i.e. straight line extrapolation.
        '''
        dz = z - DIPOLE_EDGE 
    elif z <= DIPOLE_EDGE_LOW :
        
        '''
            If the extrapolation point is upstream of the target, do something
            similar as above.
        '''
        position = extrapolate_helix_to_x_plane(track, DIPOLE_EDGE_LOW)
        dz = z - position[0]

    else: 

        '''
            If the extrapolation point is inside of the field region, 
            analytically extrapolate the helix and return the position.
        '''

        position = extrapolate_helix_to_x_plane(track, z)
        lab_position = [0]*3
        lab_position[0] = position[1]
        lab_position[1] = position[2]
        lab_position[2] = position[0]
        return lab_position

    # Calculate the azimuthal angle at the track position
    phi = get_phi(track, position)

    # Calculate the distance to the extrapolated point
    r = dz/(get_sin_theta(track)*math.cos(phi))

    # Get the delta x and y values at the point of extrapolation
    dx = r*get_sin_theta(track)*math.sin(phi)
    dy = r*get_cos_theta(track)

    # Calculate the position of the track at the extrapolation point
    x = position[1] + dx
    y = position[2] + dy

    return x, y, z
