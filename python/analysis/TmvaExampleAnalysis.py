
import ROOT
import numpy as n

signal_file = ROOT.TFile("signal.root", "RECREATE")
signal_tree = ROOT.TTree("signal", "Signal Tree")

bkg_file = ROOT.TFile("background.root", "RECREATE")
bkg_tree = ROOT.TTree("background", "Background Tree")


signal_x = n.zeros(1, dtype=float)
signal_y = n.zeros(1, dtype=float)

bkg_x = n.zeros(1, dtype=float)
bkg_y = n.zeros(1, dtype=float)

signal_tree.Branch('x', signal_x, 'x/D')
signal_tree.Branch('y', signal_y, 'y/D')

bkg_tree.Branch('x', bkg_x, 'x/D')
bkg_tree.Branch('y', bkg_y, 'y/D')

for signal_value in xrange(100000) : 
    signal_x[0] = ROOT.gRandom.Gaus(1,1)
    signal_y[0] = ROOT.gRandom.Gaus(1,1)
    signal_tree.Fill()
    
    bkg_x[0] = ROOT.gRandom.Gaus(-1, 1)
    bkg_y[0] = ROOT.gRandom.Gaus(-1, 1)
    bkg_tree.Fill()


signal_tree.SetMarkerColor(ROOT.kRed)
signal_tree.Draw("signal_y:signal_x")

bkg_tree.SetMarkerColor(ROOT.kBlue)
bkg_tree.Draw("bkg_y:bkg_x", "", "same")

signal_file.Write()

bkg_file.Write()

ROOT.TMVA.Tools.Instance()

tmva_file = ROOT.TFile("tmva_example.root", "RECREATE")

factory = ROOT.TMVA.Factory("TMVAExampleAnalysis", tmva_file, ":".join(["!V",
                                                                        "!Silent",
                                                                        "Color",
                                                                        "DrawProgressBar", 
                                                                        "Transformations=I;D;P;G;D", 
                                                                        "AnalysisType=Classification"]))


factory.AddVariable("x", "F")
factory.AddVariable("y", "F")

factory.AddSignalTree(signal_tree)
factory.AddBackgroundTree(bkg_tree)

signal_cuts = ROOT.TCut("")
bkg_cuts = ROOT.TCut("")

factory.PrepareTrainingAndTestTree(signal_cuts, bkg_cuts, ":".join(["nTrain_Signal=0", 
                                                                    "nTrain_Background=0", 
                                                                    "SplitMode=Random", 
                                                                    "NormMode=NumEvents", 
                                                                    "!V"]))


method = factory.BookMethod(ROOT.TMVA.Types.kBDT, "BDT", ":".join(["!H",
                                                                   "!V", 
                                                                   "NTrees=850", 
                                                                   "MinNodeSize=0.15%", 
                                                                   "MaxDepth=3", 
                                                                   "BoostType=AdaBoost", 
                                                                   "SeparationType=GiniIndex", 
                                                                   "nCuts=20", 
                                                                   "PruneMethod=NoPruning", 
                                                                  ]))


factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

