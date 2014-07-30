#
#
#
#
#

import math

class TrackExtrapolator:

    d0_index = 0
    phi0_index = 1
    omega_index = 2
    z0_index = 3
    tanlambda_index = 4

    def __init__(self, track_parameters):
        
        self.track_parameters = track_parameters

    def get_d0(self):
        
        return self.track_parameters[TrackExtrapolator.d0_index]

    def get_phi0(self):
        
        return self.track_parameters[self.phi0_index]

    def get_omega(self):
        
        return self.track_parameters[self.omega_index]

    def get_Z0(self):
        
        return self.track_parameters[self.z0_index]

    def get_tanLambda(self):
        
        return self.track_parameters[self.tanlambda_index]

    def get_x0(self): 
        
        return -self.get_d0()*math.sin(self.get_phi0())

    def get_y0(self):
        
        return self.get_d0()*math.cos(self.get_phi0())

    def get_R(self):
        
        return 1/self.get_omega()

    def get_xc(self):
        
        return (self.getR() - self.get_d0())*math.sin(self.get_phi0())

    def get_yc(self):
        
        return -(self.getR() - self.get_d0())*math.cos(self.get_phi0())

    def get_path_to_x_plane(self, x):
        
        r = self.get_R()
        y = self.get_yc() + math.copysign(1, r)*math.sqrt(r*r - math.pow(x - self.get_xc(), 2))

        return self.get_path_length(self.get_x0(), self.get_y0(), x, y)
        

    def get_path_length(self, x1, y1, x2, y2):
        
        phi1 = math.atan2(y1 - self.get_yc(), x1 - self.get_xc())
        phi2 = math.atan2(y2 - self.get_yc(), x2 - self.get_xc())
        dphi = phi2 - phi1

        if dphi > math.pi: dphi -= 2.*math.pi
        if dphi < -math.pi: dphi += 2.*math.pi

        return -self.get_R()*dphi

    def get_point_on_helix(self, path_length):
        
        r = self.get_R()
        phi = self.get_phi0() - path_length/r

        x = self.get_xc() - r*math.sin(phi)
        y = self.get_yc() - r*math.cos(phi)
        z = self.get_Z0() + path_length*self.get_tanLambda()

        return x, y, z

    def extrapolate_helix_to_x_plane(self, x): 
        
        path_length = self.get_path_to_x_plane(x)
        return self.get_point_on_helix(path_length)
