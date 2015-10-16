'''
    A module which provides utilties used in matching an HPS SVT track to an
    Ecal cluster.

    Author : Omar Moreno <omoreno1@ucsc.edu>
    Institution : Santa Cruz Institute for Particle Physics
                  University of California, Santa Cruz
    Data : October 15, 2015
'''

import os
import math
import TrackExtrapolator as te
import numpy as np
import matplotlib.pyplot as plt
import ROOT as r

class TrackClusterMatcher :

    def __init__(self) :

        self.top_cluster_track_match_delta_x_low = -6.10094
        self.top_cluster_track_match_delta_x_high = 12.93466
        self.bottom_cluster_track_match_delta_x_low = -8.024475
        self.bottom_cluster_track_match_delta_x_high = 10.835475
        self.bottom_cluster_track_match_delta_y_high = 7.3967
        self.bottom_cluster_track_match_delta_y_low = -8.3093
        self.top_cluster_track_match_delta_y_high = 11.494
        self.top_cluster_track_match_delta_y_low = -6.07562

        self.clusters = dict()
        self.tracks = dict()

        self.histograms = dict()
        self.histograms["cluster x - track x @ Ecal - top - all"] = r.TH1F("cluster x - track x @ Ecal - top - all",
                                                                      "cluster x - track x @ Ecal - top - all",
                                                                      200, -200, 200)
        self.histograms["cluster x - track x @ Ecal - bottom - all"] = r.TH1F("cluster x - track x @ Ecal - bottom - all",
                                                                         "cluster x - track x @ Ecal - bottom - all",
                                                                         200, -200, 200)

    def find_all_matches(self, event) :

        self.clusters.clear()
        self.tracks .clear()

        for track_n in xrange(0, event.getNumberOfTracks()) :
            track = event.getTrack(track_n)

            for cluster_n in xrange(0, event.getNumberOfEcalClusters()) :

                cluster = event.getEcalCluster(cluster_n)

                r = 0
                if self.is_match(cluster, track, r) :

                    self.clusters[cluster] = track
                    self.tracks[track] = cluster

    def is_match(self, cluster, track, r) :

        '''
           Check that the track and cluster are in the same detector volume.
           If not, there is no way they can match.
        '''
        if ((track.isTopTrack() and cluster.getPosition()[1] < 0) 
            or (track.isBottomTrack() and cluster.getPosition()[1] > 0)) : return False

        # Get the position of the cluster at the calorimeter
        cluster_pos = cluster.getPosition()

        # Extrapolate the track to the Ecal cluster position
        track_pos_at_ecal = te.extrapolate_track(track, cluster_pos[2])

        # Get the residual between the cluster and track position
        delta_x = cluster_pos[0] - track_pos_at_ecal[0]
        delta_y = cluster_pos[1] - track_pos_at_ecal[1]
        r = math.sqrt(delta_x*delta_x + delta_y*delta_y)

        if track.isTopTrack() :
            self.histograms["cluster x - track x @ Ecal - top - all"].Fill(delta_x)
        else:
            self.histograms["cluster x - track x @ Ecal - bottom - all"].Fill(delta_x)

        # Check that the track-cluster position residual is reasonable.
        if ((track.isTopTrack() and (delta_x > self.top_cluster_track_match_delta_x_high or
                                     delta_x < self.top_cluster_track_match_delta_x_low)) or
            (track.isBottomTrack() and (delta_x > self.bottom_cluster_track_match_delta_x_high or
                                        delta_x < self.bottom_cluster_track_match_delta_x_low ))) : return False;

        return True

    def save_histograms(self) :
    
        file = r.TFile("track_cluster_matching_plots.root", "RECREATE")
        self.histograms["cluster x - track x @ Ecal - top - all"].Write()
        self.histograms["cluster x - track x @ Ecal - bottom - all"].Write()


