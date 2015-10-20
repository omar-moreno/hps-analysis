import os
import ROOT as r
import numpy as n
import FlatTupleMaker as ft
import TrackClusterMatcher as tcm
import AnalysisUtils as au

class MollerBkgPreprocessor : 

    def __init__(self): 
        
        # Load the HpsEvent library.
        if os.getenv('HPS_DST_PATH') is None:
            print "[ MollerPreprocessor ]: Error! Environmental variable HPS_DST_HOME is not set."
            sys.exit(2)
        hps_dst_path = os.environ['HPS_DST_PATH']
        hps_dst_path += "/build/lib/libHpsEvent.so"
        r.gSystem.Load(hps_dst_path)
    
    def preprocess(self, dst_file):
        
        root_file = r.TFile(dst_file)
        tree = root_file.Get("HPS_Event")

        from ROOT import HpsEvent

        hps_event = HpsEvent()
        b_hps_event = tree.GetBranch("Event")
        b_hps_event.SetAddress(r.AddressOf(hps_event))

        # Create the file name for the preprocessed ROOT file
        output_file_name = root_file.GetName()[str(root_file.GetName()).rindex("/") + 1:-5]
        output_file_name += "_preprocessed_bkg.root"
        
        # Add variables to the ntuple
        ft_maker = ft.FlatTupleMaker(output_file_name)
        ft_maker.add_variable("cluster_pair_energy_high")
        ft_maker.add_variable("cluster_pair_energy_low")
        ft_maker.add_variable("cluster_pair_energy_high_x")
        ft_maker.add_variable("cluster_pair_energy_low_x")

        matcher = tcm.TrackClusterMatcher()

        for entry in xrange(0, tree.GetEntries()):
            
            tree.GetEntry(entry)

            # Loop through all clusters in the event and find a 'good' pair.  
            # For now, a 'good' pair is defined as two clusters whose cluster
            # time difference is less than 1.7 ns and greater than -1.6 ns
            cluster_pair = au.get_good_cluster_pair(hps_event)
            if len(cluster_pair) != 2 : continue
            
            matcher.find_all_matches(hps_event)
            if (matcher.get_track(cluster_pair[0]) is None) or  (matcher.get_track(cluster_pair[1]) is None) : continue

            ft_maker.set_variable_value("cluster_pair_energy_high", cluster_pair[0].getEnergy())
            ft_maker.set_variable_value("cluster_pair_energy_low", cluster_pair[1].getEnergy())
            ft_maker.set_variable_value("cluster_pair_energy_high_x", cluster_pair[0].getPosition()[0])
            ft_maker.set_variable_value("cluster_pair_energy_low_x", cluster_pair[1].getPosition()[0])

            ft_maker.fill()
        
        matcher.save_histograms()
        ft_maker.close()

