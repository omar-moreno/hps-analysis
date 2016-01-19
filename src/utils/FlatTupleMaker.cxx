/*
 * @file FlatTupleMaker.cxx
 * @author Omar Moreno
 * @date January 18, 2016
 * @brief 
 *
 */

#include <FlatTupleMaker.h>

FlatTupleMaker::FlatTupleMaker(std::string file_name, std::string tree_name) { 

    file = new TFile(file_name.c_str(), "RECREATE");

    tree = new TTree(tree_name.c_str(), tree_name.c_str());   
    
}

FlatTupleMaker::~FlatTupleMaker() { 
    delete file; 
    delete tree; 
}

void FlatTupleMaker::addVariable(std::string variable_name) { 
    
    variables[variable_name] = 0; 
    tree->Branch(variable_name.c_str(), &variables[variable_name], (variable_name + "/D").c_str()); 
}

bool FlatTupleMaker::hasVariable(std::string variable_name) { 
    
    auto search = variables.find(variable_name); 
    if (search != variables.end()) return true; 

    return false; 
}

void FlatTupleMaker::close() { 
    file->Write();
    file->Close(); 
}

