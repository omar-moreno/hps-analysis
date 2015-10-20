import os
import ROOT as r
import numpy as np
import FlatTupleMaker as ft
import TrackClusterMatcher as tcm
import AnalysisUtils as au

class MollerPreprocessor : 

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
        output_file_name += "_preprocessed_signal.root"
        
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

            ft_maker.set_variable_value("cluster_pair_energy_high", cluster_pair[0].getEnergy())
            ft_maker.set_variable_value("cluster_pair_energy_low", cluster_pair[1].getEnergy())
            ft_maker.set_variable_value("cluster_pair_energy_high_x", cluster_pair[0].getPosition()[0])
            ft_maker.set_variable_value("cluster_pair_energy_low_x", cluster_pair[1].getPosition()[0])

            ft_maker.fill()
        
        ft_maker.close()

    def is_good_event(self, event) :

        if not event.isSingle1Trigger() : return False

        if not event.isSvtBiasOn() : return False

        if not event.isSvtClosed() : return False

        if event.hasSvtBurstModeNoise() : return False

        if event.hasSvtEventHeaderErrors() : return False

        return True

    def get_good_cluster_pair(self, event) :

        cluster_pair = ()

        for first_cluster_n in xrange(event.getNumberOfEcalClusters()) :
        
            first_cluster = event.getEcalCluster(first_cluster_n)

            for second_cluster_n in xrange(event.getNumberOfEcalClusters()) :

                second_cluster = event.getEcalCluster(second_cluster_n)

                if first_cluster.getPosition()[1]*second_cluster.getPosition()[1] > 0 : continue

                cluster_pair_dt = first_cluster.getClusterTime() - second_cluster.getClusterTime()
                
                if cluster_pair_dt < -1.6 or cluster_pair_dt > 1.7 : continue

                cluster_pair = (first_cluster, second_cluster)

        return cluster_pair
