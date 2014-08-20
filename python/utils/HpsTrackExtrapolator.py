#
#
#
#
#

import math

from TrackExtrapolator import TrackExtrapolator

class HpsTrackExtrapolator(TrackExtrapolator):

	dipole_edge = 997.2
	dipole_edge_low = 0.

	def get_phi(self, position):

		x = math.sin(self.get_phi0()) - self.get_omega()*(position[0] - self.get_x0())
		y = math.cos(self.get_phi0()) + self.get_omega()*(position[1] - self.get_y0())

		return math.atan2(x, y)

	def get_sin_theta(self):
		return 1/math.sqrt(1 + math.pow(self.get_tanLambda(),2))

	def get_cos_theta(self):
		return self.get_tanLambda()/math.sqrt(1 + math.pow(self.get_tanLambda(),2))

	def extrapolate_track(self, z):
		
		if z >= HpsTrackExtrapolator.dipole_edge:

			position = self.extrapolate_helix_to_x_plane(HpsTrackExtrapolator.dipole_edge)
			dz = z - HpsTrackExtrapolator.dipole_edge
		elif z <= HpsTrackExtrapolator.dipole_edge_low:

			position = self.extrapolate_helix_to_x_plane(HpsTrackExtrapolator.dipole_edge_low)
			dz = z - position[0]

		else: 

			position = self.extrapolate_helix_to_x_plane(z)
			lab_position = [0]*3
			lab_position[0] = position[1]
			lab_position[1] = position[2]
			lab_position[2] = position[0]
			return lab_position

		phi = self.get_phi(position)

		r = dz/(self.get_sin_theta()*math.cos(phi))
		dx = r*self.get_sin_theta()*math.sin(phi)
		dy = r*self.get_cos_theta()

		x = position[1] + dx
		y = position[2] + dy

		return x, y, z

