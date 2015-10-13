
import ROOT as r
import numpy as n

class FlatTupleMaker :

    def __init__(self, file_name) :
        
        # Open a ROOT file to write the flat ntuple to
        self.root_file = r.TFile(file_name, "RECREATE")
        
        # Create a ROOT TTree used to store the tuples
        self.tree = r.TTree("preprocessed", "Preprocessed Data")

        # Dictionary containing the tuple names and values
        self.variables = dict()

    def add_variable(self, variable_name) : 

        self.variables[variable_name] = n.zeros(1, dtype=float)
        self.tree.Branch(variable_name, self.variables[variable_name], variable_name + "/D")
        
    def set_variable_value(self, variable_name, value) :

        self.variables[variable_name][0] = value

    def fill(self) : 

        self.tree.Fill()

    def close(self):

        self.root_file.Write()
        self.root_file.Close()


