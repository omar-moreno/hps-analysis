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

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <ComparePlots.h>

using namespace std;

int main(int argc, char **argv) {
   
    
    string root_file_list_name;
    bool save_to_pdf = false;

    // Parse all the command line arguments. If there are no valid command line
    // arguments passed, print the usage and exit.
    static struct option long_options[] = {
        {"root_list", required_argument, 0, 'l' },
        {"pdf", no_argument, 0, 'p' },
        { 0, 0, 0, 0 }
    };

    int option_index = 0;
    int option_char; 
    while ((option_char = getopt_long(argc, argv, "l:p", long_options, &option_index)) != -1) {
        switch (option_char) {
            case 'l': 
                root_file_list_name = optarg;
                break;
            case 'p':
                save_to_pdf = true;
                break;
        }
    }

    if (root_file_list_name.empty()) { 
        cerr << "Please specify a list of ROOT files to process." << endl;
        return EXIT_FAILURE;
    }

    // Create the list of files to process
    list< TFile* > root_files;
    string root_file; 
    ifstream root_file_list(root_file_list_name.c_str(), ifstream::in);
    if (!root_file_list.is_open()) { 
        cerr << "Failed to open file " << root_file_list_name << endl;
        return EXIT_FAILURE;
    }

    while (root_file_list >> root_file) { 
        root_files.push_back(new TFile(root_file.c_str())); 
        if (root_files.back()->IsZombie()) { 
            cout << "Failed to open root file " << root_file << endl;
            return EXIT_FAILURE;
        }
    } 
    root_file_list.close();

    ComparePlots comparator;
    comparator.parseFiles(root_files);
    comparator.overlayPlots(); 

    if (save_to_pdf) comparator.saveToPdf("test.pdf");

    return EXIT_SUCCESS;
}


