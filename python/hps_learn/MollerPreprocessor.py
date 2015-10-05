import os
import ROOT as root
import numpy as n

class MollerPreprocessor : 

    def __init__(self): 
        # Load the HpsEvent library.
        if os.getenv('HPS_DST_PATH') is None:
            print "[ MollerPreprocessor ]: Error! Environmental variable HPS_DST_HOME is not set."
            sys.exit(2)
        hps_dst_path = os.environ['HPS_DST_PATH']
        hps_dst_path += "/build/lib/libHpsEvent.so"
        root.gSystem.Load(hps_dst_path)
        

    def preprocess_signal(self, dst_file):
      
        tree = dst_file.Get("HPS_Event")

        from ROOT import HpsEvent

        hps_event = HpsEvent()
        b_hps_event = tree.GetBranch("Event")
        b_hps_event.SetAddress(root.AddressOf(hps_event))

        output_file_name = dst_file.GetName()[str(dst_file.GetName()).rindex("/") + 1:-5]
        output_file_name += "_preprocessed_signal.root"

        output_file = root.TFile(output_file_name, "RECREATE")
        output_tree = root.TTree("Signal", "Signal Tree")

        first_cluster_energy = n.zeros(1, dtype=float)
        first_cluster_x_position = n.zeros(1, dtype=float)
        second_cluster_energy = n.zeros(1, dtype=float)
        second_cluster_x_position = n.zeros(1, dtype=float)

        output_tree.Branch('first_cluster_energy', first_cluster_energy, 'first_cluster_energy/D')
        output_tree.Branch('first_cluster_x_position', first_cluster_x_position, 'first_cluster_x_position/D')
        output_tree.Branch('second_cluster_energy', second_cluster_energy, 'second_cluster_energy/D')
        output_tree.Branch('second_cluster_x_position', second_cluster_x_position, 'second_cluster_x_position/D')

        for entry in xrange(0, tree.GetEntries()):

            tree.GetEntry(entry)

            if hps_event.getNumberOfEcalClusters() != 2: continue
            
            first_cluster_energy[0] = hps_event.getEcalCluster(0).getEnergy()
            first_cluster_x_position[0] = hps_event.getEcalCluster(0).getPosition()[0]
            second_cluster_energy[0] = hps_event.getEcalCluster(1).getEnergy()
            second_cluster_x_position[0] = hps_event.getEcalCluster(1).getPosition()[0]
           
            output_tree.Fill()
            
        output_file.Write()
        return output_file_name

    def preprocess_background(self, dst_file): 
        
        tree = dst_file.Get("HPS_Event")

        from ROOT import HpsEvent

        hps_event = HpsEvent()
        b_hps_event = tree.GetBranch("Event")
        b_hps_event.SetAddress(root.AddressOf(hps_event))

        output_file_name = dst_file.GetName()[str(dst_file.GetName()).rindex("/") + 1:-5]
        output_file_name += "_preprocessed_bkg.root"

        output_file = root.TFile(output_file_name, "RECREATE")
        output_tree = root.TTree("Background", "Background Tree")

        first_cluster_energy = n.zeros(1, dtype=float)
        first_cluster_x_position = n.zeros(1, dtype=float)
        second_cluster_energy = n.zeros(1, dtype=float)
        second_cluster_x_position = n.zeros(1, dtype=float)

        output_tree.Branch('first_cluster_energy', first_cluster_energy, 'first_cluster_energy/D')
        output_tree.Branch('first_cluster_x_position', first_cluster_x_position, 'first_cluster_x_position/D')
        output_tree.Branch('second_cluster_energy', second_cluster_energy, 'second_cluster_energy/D')
        output_tree.Branch('second_cluster_x_position', second_cluster_x_position, 'second_cluster_x_position/D')

        for entry in xrange(0, tree.GetEntries()):

            tree.GetEntry(entry)

            if not hps_event.isSingle1Trigger() : continue

            if hps_event.getNumberOfEcalClusters() != 2: continue

            delta_t = abs(hps_event.getEcalCluster(0).getClusterTime() - hps_event.getEcalCluster(1).getClusterTime())
            if delta_t > 3 : continue

            first_cluster_energy[0] = hps_event.getEcalCluster(0).getEnergy()
            first_cluster_x_position[0] = hps_event.getEcalCluster(0).getPosition()[0]
            second_cluster_energy[0] = hps_event.getEcalCluster(1).getEnergy()
            second_cluster_x_position[0] = hps_event.getEcalCluster(1).getPosition()[0]
           
            output_tree.Fill()
            
        output_file.Write()
        return output_file_name


