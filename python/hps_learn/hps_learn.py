#!/usr/bin/python

###############
#   Imports   #
###############

import argparse
import sys
import ROOT as root
from ConfigReader import ConfigReader

############
#   Main   #
############

def main():

    # Parse all command line arguments using the argparse module
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-c", "--config", help="Configuration file")
    args = parser.parse_args()

    if args.config is None: 
        print '[ hps-learn ]: A configuration file needs to be specified!'
        sys.exit(2)

    # Parse the configuration
    config_reader = ConfigReader(args.config)

    # If signal files were specified, open them and get the specified ROOT 
    # TTree objects.  Otherwise, warn the user that no signal files were
    signal_files = config_reader.get_signal_files()
    signal_root_files = list()
    signal_root_trees = list()
    if signal_files :
        for signal_file, tree_name in signal_files.iteritems():
            print "[ hps-learn ]: Loading signal file " + signal_file + " with TTree name " + tree_name
            signal_root_files.append(root.TFile(signal_file))
            signal_root_trees.append(signal_root_files[-1].Get(tree_name))
    else :
        print "[ hps-learn ]: At least a single signal file needs to be specified."
        sys.exit(2)
   
    # If background files were specified, open them and get the specified ROOT 
    # TTree objects. Otherwise, warn the user that no background files were 
    # specified and exit.
    bkg_files = config_reader.get_background_files()
    bkg_root_files = list()
    bkg_root_trees = list()
    if bkg_files:
        for bkg_file, tree_name in bkg_files.iteritems(): 
            print "[ hps-learn ]: Loading background file " + bkg_file + " with TTree name " + tree_name
            bkg_root_files.append(root.TFile(bkg_file))
            bkg_root_trees.append(bkg_root_files[-1].Get(tree_name))
    else : 
        print "[ hps-learn ]: At least a single background file needs to be specified."
        sys.exit(2)

    root.TMVA.Tools.Instance()

    tmva_file = root.TFile("tmva_output.root", "RECREATE")

    factory = root.TMVA.Factory("TMVA_Analysis", tmva_file, ":".join(["!V",
                                                                "!Silent",
                                                                "Color",
                                                                "DrawProgressBar", 
                                                                "Transformations=I;D;P;G;D", 
                                                                "AnalysisType=Classification"]))

    for signal_tree in signal_root_trees:
        factory.AddSignalTree(signal_tree)

    for bkg_tree in bkg_root_trees:
        factory.AddBackgroundTree(bkg_tree)

    variables = config_reader.get_training_variables()
    for variable, variable_type in variables.iteritems():
        print "Adding training variable " + variable + " of type " + variable_type
        factory.AddVariable(variable, variable_type)

    signal_cuts = root.TCut("")
    bkg_cuts = root.TCut("")

    factory.PrepareTrainingAndTestTree(signal_cuts, bkg_cuts, ":".join(["nTrain_Signal=0", 
                                                                    "nTrain_Background=0", 
                                                                    "SplitMode=Random", 
                                                                    "NormMode=NumEvents", 
                                                                    "!V"]))

    method = factory.BookMethod(root.TMVA.Types.kMLP, "MLP_ANN", "")

    factory.TrainAllMethods()
    factory.TestAllMethods()
    factory.EvaluateAllMethods()

if __name__ == "__main__":
    main()
