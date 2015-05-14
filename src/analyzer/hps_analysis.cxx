
//-----------//
//--- C++ ---//
//-----------//
#include <cstdlib>
#include <iostream>
#include <getopt.h>
#include <string>
#include <list>

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

using namespace std; 

void printUsage(); 

int main(int argc, char **argv) {
    
    string dst_file_name; 

    // Parse the command line arguments
    static struct option long_options[] = 
    {
        {"input_file", required_argument, 0, 'i'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    int option_index = 0;
    int option_char; 
    while((option_char = getopt_long(argc, argv, "i:h", long_options, &option_index)) != -1) {

        switch(option_char){
            case 'i': 
                dst_file_name = optarg;
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

    // If a DST file name is not specified, warn the user and exit the program
    if (dst_file_name.empty()) { 
        cerr << "\n[ HPS ANALYZER ]: Please specify a DST file to process." << endl;
        return EXIT_FAILURE;
    } 

    // Open the ROOT file
 	TFile *file = new TFile(dst_file_name.c_str());

    // Get the TTree "HPS_EVENT" containing the HpsEvent branch and all
    // other collections
    TTree *tree = (TTree*) file->Get("HPS_Event");

    // Create a pointer to an HpsEvent object in order to read the TClonesArrays
    // collections
    HpsEvent *event = new HpsEvent();

    // Get the HpsEvent branch from the TTree and set the branch address to
    // the pointer created above
    TBranch *b_event = tree->GetBranch("Event");
    b_event->SetAddress(&event);

    // Container to hold all analyses
    list<HpsAnalysis*> analyses;

    analyses.push_back(new TrackAnalysis());

    // Initialize all analyses
    for (list<HpsAnalysis*>::iterator analysis = analyses.begin();
            analysis != analyses.end(); ++analysis) { 
        cout << "[ HPS ANALYZER ]: Initializing analysis: " << (*analysis)->toString() << endl;
        (*analysis)->initialize();
    }

    for (int entry = 0; entry < tree->GetEntries(); ++entry) {
   
    	// Print the event number every 500 events
    	if((entry+1)%500 == 0){
    		std::cout << "Event: " << entry+1 << endl;
    	}
        
        // Read the ith entry from the tree.  This "fills" HpsEvent and allows
        // access to all collections
        tree->GetEntry(entry);

        for (list<HpsAnalysis*>::iterator analysis = analyses.begin();
                analysis != analyses.end(); ++analysis) { 
            (*analysis)->processEvent(event);
        }
    }

    for (list<HpsAnalysis*>::iterator analysis = analyses.begin();
            analysis != analyses.end(); ++analysis) { 
        (*analysis)->finalize();
        delete *analysis;
    }
    analyses.clear();

    return EXIT_SUCCESS;  
}

void printUsage(){}
