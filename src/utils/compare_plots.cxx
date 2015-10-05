/**
 *  @file ComparePlots.h
 *  @brief 
 *  @author Omar Moreno <omoreno1@ucsc.edu>
 *  @date April 29, 2015
 *
 */

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <getopt.h>
#include <iostream>
#include <cstdlib>
#include <list>
#include <fstream>

//--------------//
//   yaml-cpp   //
//--------------//
#include <yaml-cpp/yaml.h>

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <ComparePlots.h>

using namespace std;

int main(int argc, char **argv) {
   
    
    string config_file_name; 
    string style = "";
    bool save_to_pdf = false;

    // Parse all the command line arguments. If there are no valid command line
    // arguments passed, print the usage and exit.
    static struct option long_options[] = {
        {"config_file", required_argument, 0, 'c' }, 
        {"style", required_argument, 0, 's'},
        {"pdf", no_argument, 0, 'p' },
        { 0, 0, 0, 0 }
    };

    int option_index = 0;
    int option_char; 
    while ((option_char = getopt_long(argc, argv, "c:s:p", long_options, &option_index)) != -1) {
        switch (option_char) {
            case 'c': 
                config_file_name = optarg; 
                break;
            case 's':
                style = optarg; 
                break;
            case 'p':
                save_to_pdf = true;
                break;
            default: 
                break;
        }
    }

    // Require a configuration
    if (config_file_name.empty()) { 
        cerr << "Please specify a configuration file." << endl;
        return EXIT_FAILURE; 
    }

    // TODO: This all needs to move to it's own config parse method
    YAML::Node config = YAML::LoadFile(config_file_name);
    if (config.IsNull()) { 
        cerr << "Unable to open configuration file " << config_file_name << std::endl;
        return EXIT_FAILURE;
    } 
    
    YAML::Node root_file_names = config["Files"];
    if (root_file_names.IsNull()) { 
        cerr << "Please specify a list of ROOT files to process." << endl;
        return EXIT_FAILURE;    
    }

    list< TFile* > root_files;
    for (auto root_file_name : root_file_names) { 
        root_files.push_back(new TFile(root_file_name.as<string>().c_str())); 
        if (root_files.back()->IsZombie()) { 
            cout << "Failed to open ROOT file " << root_file_name << std::endl;
            return EXIT_FAILURE; 
        } 
    }

    ComparePlots* comparator = new ComparePlots();
    if (!style.empty()) comparator->setStyle(style);

    YAML::Node histogram_names = config["Histograms"];
    if (!histogram_names.IsNull()) { 
        comparator->parseFiles(root_files, histogram_names[0].as<string>().c_str());
    } else comparator->parseFiles(root_files);
    comparator->overlayPlots(); 

    return EXIT_SUCCESS;
}


