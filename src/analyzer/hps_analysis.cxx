
//--- C++ ---//
//-----------//
#include <cstdlib>
#include <iostream>
#include <getopt.h>
#include <string>

//--- libxml2 ---//
//---------------//
#include <libxml/parser.h>

//--- ROOT ---//
//------------//

using namespace std; 

void printUsage(); 

int main(int argc, char **argv)
{

    
    string config_file_name; 

    // Parse the command line arguments
    static struct option long_options[] = 
    {
        {"config_file", required_argument, 0, 'i'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    int option_index = 0;
    int option_char; 
    while((option_char = getopt_long(argc, argv, "i:h", long_options, &option_index)) != -1){

        switch(option_char){
            case 'i': 
                config_file_name = optarg;
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
    if(config_file_name.empty()){
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

   return EXIT_SUCCESS;  
}

void printUsage(){}
