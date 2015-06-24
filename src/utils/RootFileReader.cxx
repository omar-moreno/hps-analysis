/**
 *  @file RootFileReader.cxx
 *  @brief Reader used to parse ROOT files. 
 *  @author Omar Moreno <omoreno1@ucsc.edu>
 *  @date June 22, 2015
 */

#include <RootFileReader.h>

RootFileReader::RootFileReader()
    : histogram_substring("")
{
}

RootFileReader::~RootFileReader() { 
}

void RootFileReader::parseFile(TFile* root_file) { 

    std::cout << "[ RootFileReader ]: Processing file: " << root_file->GetName() << std::endl;
    
    this->parseFile(root_file->GetListOfKeys());
}


void RootFileReader::parseFile(TList* keys) { 

    // Instantiate an iterator that will be used to loop through all of the 
    // histograms in the file.
    TIter next(keys);

    // Iterate through the histograms (keys)    
    while (TKey *key = (TKey*) next()) { 

        // If the key points to a directory, retrieve the keys inside of the 
        // directory and load them to the corresponding map
        if (key->IsFolder()) this->parseFile(((TDirectory*) key->ReadObj())->GetListOfKeys());
        
        // 
        if (!histogram_substring.empty() 
                && std::string(key->GetName()).find(histogram_substring) == std::string::npos)
            continue;

        //std::cout << "Object type: " << key->ReadObj()->ClassName() << std::endl;
        if (std::string(key->ReadObj()->ClassName()).find("1") != std::string::npos) {
            histogram1D_map[key->GetName()].push_back((TH1*) key->ReadObj()); 
        } else if (std::string(key->ReadObj()->ClassName()).find("2") != std::string::npos) {
            histogram2D_map[key->GetName()].push_back((TH1*) key->ReadObj()); 
        } else if (std::string(key->ReadObj()->ClassName()).find("Graph") != std::string::npos) { 
            graph_map[key->GetName()].push_back((TGraph*) key->ReadObj());
        }
        std::cout << "[ RootFileReader ]: Adding " << key->GetName() << std::endl;
    } 
}

void RootFileReader::parseFile(TFile* root_file, std::string histogram_substring) { 
    
    histogram1D_map.clear();
    histogram2D_map.clear();

    this->setHistogramName(histogram_substring);
    
    this->parseFile(root_file);
}

void RootFileReader::parseFiles(std::list<TFile*> root_files) { 

    // Loop over all of the ROOT files and create the plot maps
    std::list<TFile*>::iterator root_files_it = root_files.begin();
    for (root_files_it; root_files_it != root_files.end(); ++root_files_it) { 
        std::cout << "[ ComparePlots ]: Processing file " << (*root_files_it)->GetName() << std::endl;        
        this->parseFile(*root_files_it);    
    }
} 

void RootFileReader::parseFiles(std::list<TFile*> root_files, std::string histogram_substring) { 

    // Loop over all of the ROOT files and create the plot maps
    std::list<TFile*>::iterator root_files_it = root_files.begin();
    for (root_files_it; root_files_it != root_files.end(); ++root_files_it) { 
        std::cout << "[ ComparePlots ]: Processing file " << (*root_files_it)->GetName() << std::endl;        
        this->parseFile(*root_files_it, histogram_substring);    
    }
} 

std::vector<TH1*> RootFileReader::getMatching1DHistograms(std::string histogram_name) { 

    std::vector<TH1*> histogram_collection;

    // Iterate through the collection of 1D histograms and add those matching
    // histogram_name to the list
    std::map<std::string, std::vector<TH1*>>::iterator histogram1D_it = histogram1D_map.begin();
    for (histogram1D_it; histogram1D_it != histogram1D_map.end(); histogram1D_it++) {
        
       if (histogram1D_it->first.find(histogram_name) != std::string::npos) {
            histogram_collection.insert(histogram_collection.begin(), 
                    histogram1D_it->second.begin(), histogram1D_it->second.end());
       }
    }

   return histogram_collection; 
}
