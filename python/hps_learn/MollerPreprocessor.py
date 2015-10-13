import os
import ROOT as r
import numpy as n
import FlatTupleMaker as ft
import EcalUtilsModule as em
import ctypes

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

            ft_maker.set_variable_value("cluster_pair_energy_high", 10)
            ft_maker.set_variable_value("cluster_pair_energy_low", 10)

            ft_maker.fill()

            '''

            if hps_event.getNumberOfEcalClusters() != 2: continue
            
            first_cluster_energy[0] = hps_event.getEcalCluster(0).getEnergy()
            first_cluster_x_position[0] = hps_event.getEcalCluster(0).getPosition()[0]
            second_cluster_energy[0] = hps_event.getEcalCluster(1).getEnergy()
            second_cluster_x_position[0] = hps_event.getEcalCluster(1).getPosition()[0]
           
            output_tree.Fill()
            '''
        ft_maker.close()

    def get_good_cluster_pair(self, event) :

        for first_cluster_n in xrange(event.getEcalCluster()) :
        
            first_ckuster = event.getEcalCluster(first_cluster_n)

            for second_cluster_n in xrange(event.getEcalCluster()) :

                second_cluster = event.getEcalCluster(second_cluster_n)
