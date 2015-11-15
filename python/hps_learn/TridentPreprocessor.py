import os
import sys
import math
import ROOT as r
import numpy as np
import FlatTupleMaker as ft
import TrackClusterMatcher as tcm
import TrackExtrapolator as te
import AnalysisUtils as au

class TridentPreprocessor(object) : 

    def __init__(self): 
        # Load the HpsEvent library.
        if os.getenv('HPS_DST_PATH') is None:
            print "[ MollerPreprocessor ]: Error! Environmental variable HPS_DST_HOME is not set."
            sys.exit(2)
        hps_dst_path = os.environ['HPS_DST_PATH']
        hps_dst_path += "/build/lib/libHpsEvent.so"
        r.gSystem.Load(hps_dst_path)

        self.output_file_name = "bh_preprocessed_file_test.root"

    def preprocess(self, dst_file):
      
        chain = r.TChain("HPS_Event")
        chain.Add(dst_file)
        tree = chain.GetTree()

        from ROOT import HpsEvent

        hps_event = HpsEvent()
        branch = chain.SetBranchAddress("Event", r.AddressOf(hps_event))

        # Add variables to the ntuple
        ft_maker = ft.FlatTupleMaker(self.output_file_name)
        ft_maker.add_variable("cluster_0_energy")
        ft_maker.add_variable("cluster_1_energy")
        ft_maker.add_variable("cluster_0_x")
        ft_maker.add_variable("cluster_1_x")
        ft_maker.add_variable("cluster_0_y")
        ft_maker.add_variable("cluster_1_y")
        ft_maker.add_variable("track_0_p")
        ft_maker.add_variable("track_0_px")
        ft_maker.add_variable("track_0_py")
        ft_maker.add_variable("track_0_pz")
        ft_maker.add_variable("track_0_theta")
        ft_maker.add_variable("track_0_phi0")
        ft_maker.add_variable("track_0_omega")
        ft_maker.add_variable("track_0_d0")
        ft_maker.add_variable("track_0_z0")
        ft_maker.add_variable("track_1_p")
        ft_maker.add_variable("track_1_px")
        ft_maker.add_variable("track_1_py")
        ft_maker.add_variable("track_1_pz")
        ft_maker.add_variable("track_1_theta")
        ft_maker.add_variable("track_1_phi0")
        ft_maker.add_variable("track_1_omega")
        ft_maker.add_variable("track_1_d0")
        ft_maker.add_variable("track_1_z0")
        ft_maker.add_variable("track_pair_p_sum")

        matcher = tcm.TrackClusterMatcher()

        for entry in xrange(0, chain.GetEntries()):
           
            if (entry+1)%1000 == 0 : print "Event " + str(entry+1)

            chain.GetEntry(entry)
        
            # Loop through all clusters in the event and find a 'good' pair.  
            # For now, a 'good' pair is defined as two clusters whose cluster
            # time difference is less than 1.6 ns and greater than -1.6 ns
            cluster_pair = au.get_good_cluster_pair(hps_event)

            if len(cluster_pair) != 2 : continue 

            matcher.find_all_matches(hps_event)
            tracks = [matcher.get_track(cluster_pair[0]), matcher.get_track(cluster_pair[1])]
            
            if (tracks[0] is None) or (tracks[1] is None) : continue

            ft_maker.set_variable_value("cluster_0_energy", cluster_pair[0].getEnergy())
            ft_maker.set_variable_value("cluster_1_energy", cluster_pair[1].getEnergy())
            ft_maker.set_variable_value("cluster_0_x", cluster_pair[0].getPosition()[0])
            ft_maker.set_variable_value("cluster_1_x", cluster_pair[1].getPosition()[0])
            ft_maker.set_variable_value("cluster_0_y", cluster_pair[0].getPosition()[1])
            ft_maker.set_variable_value("cluster_1_y", cluster_pair[1].getPosition()[1])
            ft_maker.set_variable_value("track_0_p", np.linalg.norm(np.asarray(tracks[0].getMomentum())))
            ft_maker.set_variable_value("track_0_px", tracks[0].getMomentum()[0])
            ft_maker.set_variable_value("track_0_py", tracks[0].getMomentum()[1])
            ft_maker.set_variable_value("track_0_pz", tracks[0].getMomentum()[2])
            track_0_theta = math.fabs(math.pi/2 - math.acos(te.get_cos_theta(tracks[0])))
            ft_maker.set_variable_value("track_0_theta", track_0_theta)
            ft_maker.set_variable_value("track_0_phi0", tracks[0].getPhi0())
            ft_maker.set_variable_value("track_0_omega", tracks[0].getOmega())
            ft_maker.set_variable_value("track_0_d0", tracks[0].getD0())
            ft_maker.set_variable_value("track_0_z0", tracks[0].getZ0())
            ft_maker.set_variable_value("track_1_p", np.linalg.norm(np.asarray(tracks[1].getMomentum())))
            ft_maker.set_variable_value("track_1_px", tracks[1].getMomentum()[0])
            ft_maker.set_variable_value("track_1_py", tracks[1].getMomentum()[1])
            ft_maker.set_variable_value("track_1_pz", tracks[1].getMomentum()[2])
            track_1_theta = math.fabs(math.pi/2 - math.acos(te.get_cos_theta(tracks[1])))
            ft_maker.set_variable_value("track_1_theta", track_1_theta)
            ft_maker.set_variable_value("track_1_phi0", tracks[1].getPhi0())
            ft_maker.set_variable_value("track_1_omega", tracks[1].getOmega())
            ft_maker.set_variable_value("track_1_d0", tracks[1].getD0())
            ft_maker.set_variable_value("track_1_z0", tracks[1].getZ0())
            p_sum = np.linalg.norm(np.asarray(tracks[0].getMomentum())) + np.linalg.norm(np.asarray(tracks[1].getMomentum()))
            ft_maker.set_variable_value("track_pair_p_sum", p_sum)
            
            ft_maker.fill()
        
        ft_maker.close()

    def set_output_file_name(self, output_file_name) : 
        self.output_file_name = output_file_name
