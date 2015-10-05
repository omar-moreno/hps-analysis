/**
 *  @file RootFileReader.cxx
 *  @brief Reader used to parse ROOT files. 
 *  @author Omar Moreno <omoreno1@ucsc.edu>
 *  @date June 22, 2015
 */

#include <RootFileReader.h>

RootFileReader::RootFileReader()
    : histogram_substring("") {
}

RootFileReader::~RootFileReader() { 
}

void RootFileReader::parseFile(TFile* root_file) { 


    std::cout << "[ RootFileReader ]: Processing file: " << root_file->GetName() << std::endl; 
    this->parseFile(root_file->GetListOfKeys());
}


void RootFileReader::parseFile(TList* keys) { 

    // Instantiate an iterator that will be used to loop through all of the
    // keys retrieved from the ROOT file.
    TIter next(keys);

    // Iterate through all of the keys in the file.  A key can be 
    // representitive of a directory, histogram or a graph.  
    while (TKey *key = (TKey*) next()) { 

        // If the key is associated with a directory, enter the directory, 
        // retrieve all of the keys and iterate through them. This is done so
        // the reader can handle any tree structure. 
        if (key->IsFolder()) this->parseFile(((TDirectory*) key->ReadObj())->GetListOfKeys());
        
        // If a user has specified a subset of histograms or graphs by name, 
        // check that the histogram/graph associated with the current key is 
        // named accordingly.
        if (!histogram_substring.empty() 
                && std::string(key->GetName()).find(histogram_substring) == std::string::npos) continue;

        // Sort the histograms based on type
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
   
    // Add the names of the histograms of interest. 
    this->setHistogramName(histogram_substring);
    
    // Parse the ROOT file and extract the histograms of interest
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

std::vector<TH1*> RootFileReader::getMatching2DHistograms(std::string histogram_name) { 

    std::vector<TH1*> histogram_collection;

    // Iterate through the collection of 1D histograms and add those matching
    // histogram_name to the list
    std::map<std::string, std::vector<TH1*>>::iterator histogram2D_it = histogram2D_map.begin();
    for (histogram2D_it; histogram2D_it != histogram2D_map.end(); histogram2D_it++) {
        
       if (histogram2D_it->first.find(histogram_name) != std::string::npos) {
            histogram_collection.insert(histogram_collection.begin(), 
                    histogram2D_it->second.begin(), histogram2D_it->second.end());
       }
    }

   return histogram_collection; 
}

std::vector<TGraph*> RootFileReader::getMatchingGraphs(std::string graph_name) { 

    std::vector<TGraph*> graph_collection;

    std::map<std::string, std::vector<TGraph*>>::iterator graph_it = graph_map.begin();
    for (graph_it; graph_it != graph_map.end(); graph_it++) {
        
       if (graph_it->first.find(graph_name) != std::string::npos) {
            graph_collection.insert(graph_collection.begin(), 
                    graph_it->second.begin(), graph_it->second.end());
       }
    }

   return graph_collection; 
}
