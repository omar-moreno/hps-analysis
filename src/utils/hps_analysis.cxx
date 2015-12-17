/**
 * @file hps_analsysis.cxx
 * @brief Application used to process and analyze DST events.
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date 
 */


//------------------//
//--- C++ StdLib ---//
//------------------//
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <string>
#include <list>

//---------------//
//--- libxml2 ---//
//---------------//
//#include <libxml/parser.h>

//------------//
//--- ROOT ---//
//------------//
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <TrackAnalysis.h>
#include <TagProbeAnalysis.h>
#include <TrackClusterMatchingEfficiencyAnalysis.h>
#include <SharedHitAnalysis.h>
#include <TimingAnalysis.h>
#include <MollerAnalysis.h>
#include <MuonAnalysis.h>
#include <GblTrackAnalysis.h>
#include <TridentAnalysis.h>
#include <TridentDataAnalysis.h>
#include <MollerDataAnalysis.h>

using namespace std; 

void printUsage(); 

int main(int argc, char **argv) {

    string dst_file_name;
    string file_list_name;
    int event_count = -1; 

    // Parse all the command line arguments.  If there are no valid command
    // line arguments passed, print the usage and exit the application
    static struct option long_options[] = {
        {"file_name",  required_argument, 0, 'i'},
        {"file_list",  required_argument, 0, 'l'},
        {"events",     required_argument, 0, 'n'},
        {"help",       no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int option_char; 
    while ((option_char = getopt_long(argc, argv, "i:l:n:h", long_options, &option_index)) != -1) {

        switch(option_char){
            case 'i': 
                dst_file_name = optarg;
                break;
            case 'l':
                file_list_name = optarg;
                break; 
            case 'n':
                event_count = atoi(optarg);
                break;
            case 'h':
                printUsage();
                return EXIT_SUCCESS; 
            default: 
                printUsage(); 
                return EXIT_FAILURE;
        }
    }

    // If an input file name is not specified, exit the program
    /*if(config_file_name.empty()){
      cerr << "\nPlease specify a configuration file.\n" 
      << "Use --help for usage\n" << endl;
      return EXIT_FAILURE; 
      }

      LIBXML_TEST_VERSION

    // Parse the input file and extract the configuration
    xmlDoc* doc = NULL; 

    doc = xmlReadFile(config_file_name.c_str(), NULL, 0);
    if(doc == NULL){
    cerr << "Unable to open file " << config_file_name << endl; 
    return EXIT_FAILURE; 
    } 

    xmlNode* root_element = NULL; 
    root_element = xmlDocGetRootElement(doc);
    cout << "Root element name: " << root_element->name << endl;  

    xmlFreeDoc(doc);
    xmlCleanupParser();  
    */

    // If a DST file or a list of files was not specified, warn the user 
    // and exit the application.  Also exit if both a file and a list of files
    // has been specified
    if (dst_file_name.empty() && file_list_name.empty()) { 
        cerr << "\n[ HPS ANALYZER ]: Please specify a DST file to process." << endl;
        cerr << "[ HPS ANALYZER ]: Use --help for usage.\n" << endl;
        return EXIT_FAILURE;
    } else if (!dst_file_name.empty() && !file_list_name.empty()) { 
        cerr << "\n[ HPS ANALYZER ]: Cannot specify both an DST file name and a "
             << "list of files." << endl;
        cerr << "[ HPS ANALYZER ]: Use --help for usage.\n" << endl;
        return EXIT_FAILURE;
    } 

    // Create a list of files to process
    list<string> files; 
    string file;
    if (!dst_file_name.empty()) { 
        files.push_back(dst_file_name); 
    } else if (!file_list_name.empty()) { 
        
        ifstream file_list(file_list_name.c_str(), ifstream::in);
        if (!file_list.is_open()) { 
            cerr << "\n[ HPS ANALYZER ]: Failed to open file " << file_list_name << endl;
            return EXIT_FAILURE;
        }
        
        while (file_list >> file) { 
            files.push_back(file); 
        }
        file_list.close();
    }

    // Container to hold all analyses
    list<HpsAnalysis*> analyses;

    //analyses.push_back(new TrackAnalysis());
    //analyses.push_back(new GblTrackAnalysis());
    //analyses.push_back(new TimingAnalysis());
    //analyses.push_back(new TagProbeAnalysis());
    //analyses.push_back(new TrackClusterMatchingEfficiencyAnalysis());
    //analyses.push_back(new SharedHitAnalysis());
    //analyses.push_back(new MollerAnalysis());
    //analyses.push_back(new MollerDataAnalysis());
    //analyses.push_back(new MuonAnalysis());
    analyses.push_back(new TridentAnalysis()); 
    //analyses.push_back(new TridentDataAnalysis()); 

    // Initialize all analyses
    for (list<HpsAnalysis*>::iterator analysis = analyses.begin();
            analysis != analyses.end(); ++analysis) { 
        cout << "[ HPS ANALYZER ]: Initializing analysis: " << (*analysis)->toString() << endl;
        (*analysis)->initialize();
    }

    // Create a pointer to an HpsEvent object in order to read the TClonesArrays
    // collections
    HpsEvent *event = new HpsEvent();

    // Loop over all input files and process them
    for (list<string>::iterator files_it = files.begin(); files_it != files.end(); ++files_it) { 
        
        // Open the ROOT file.  If the file can't be opened, exit the 
        // application
        TFile* file = new TFile((*files_it).c_str());

        cout << "[ HPS ANALYZER ]: Processing file: " << *files_it << endl;

        // Get the TTree "HPS_EVENT" containing the HpsEvent branch and all
        // other collections
        TTree* tree = (TTree*) file->Get("HPS_Event");
    
        // Get the HpsEvent branch from the TTree and set the branch address to
        // the pointer created above
        TBranch *b_event = tree->GetBranch("Event");
        b_event->SetAddress(&event);

        // Loop over all of the events and process them
        for (int entry = 0; entry < tree->GetEntries(); ++entry) {
        
            // Print the event number every 500 events
            if((entry+1)%10000 == 0){
                std::cout << "[ HPS ANALYZER ]: Event: " << entry+1 << endl;
            }
        
            // Read the ith entry from the tree.  This "fills" HpsEvent and allows
            // access to all collections
            tree->GetEntry(entry);

            // Loop over all analyses in the list and process the events
            for (list<HpsAnalysis*>::iterator analysis = analyses.begin();
                    analysis != analyses.end(); ++analysis) { 
                (*analysis)->processEvent(event);
            }
        
            if (entry+1 == event_count) break;
        }
        
        // Delete the file
        delete file;
    }

    // Finalize all of the analyses and free them
    for (list<HpsAnalysis*>::iterator analysis = analyses.begin();
            analysis != analyses.end(); ++analysis) { 
        (*analysis)->finalize();
        delete *analysis;
    }
    analyses.clear();

    return EXIT_SUCCESS;  
}

void printUsage() {
    cout << "\nUsage: hps_analysis [OPTIONS]" << endl;
}
