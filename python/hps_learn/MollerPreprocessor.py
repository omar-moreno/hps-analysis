import os
import ROOT as r
import numpy as n
import FlatTupleMaker as ft

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

        for entry in xrange(0, tree.GetEntries()):
            
            tree.GetEntry(entry)

            cluster_pair = self.get_good_cluster_pair(hps_event)
            print len(cluster_pair)

            ft_maker.set_variable_value("cluster_pair_energy_high", 10)
            ft_maker.set_variable_value("cluster_pair_energy_low", 10)

            ft_maker.fill()
        
        ft_maker.close()

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
