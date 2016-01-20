/**
 *
 */

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <cstdlib>
#include <iostream>
#include <getopt.h>
#include <string>

#include <TFile.h>
#include <TH1.h>
#include <RooFitResult.h>

#include <RootFileReader.h>
#include <BumpHunter.h>
#include <FlatTupleMaker.h>

using namespace std;

int main(int argc, char **argv) { 

    string file_name;
    int hist_count = 10;
    int poly_order = 3;
    double window_size = 0.020; 
    double window_start = 0.03;
    double window_end = 0.050;
    bool bkg_only = false; 

    // Parse all the command line arguments.  If there are no valid command
    // line arguments passed, print the usage and exit the application
    static struct option long_options[] = {
        {"bkg_only",   required_argument, 0, 'b'},
        {"file_name",  required_argument, 0, 'i'},
        {"number",     required_argument, 0, 'n'},
        {"order",      required_argument, 0, 'o'},
        {"window",     required_argument, 0, 'w'},
        {"start",      required_argument, 0, 's'},
        {"end",        required_argument, 0, 'e'},
        {"help",       no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };
    
    int option_index = 0;
    int option_char; 
    while ((option_char = getopt_long(argc, argv, "bi:n:o:w:s:e:h", long_options, &option_index)) != -1) {

        switch(option_char) {
            case 'b': 
                bkg_only = true;
                break; 
            case 'i': 
                file_name = optarg;
                break;
            case 'n':
                hist_count = atoi(optarg);
                break;
            case 'o': 
                poly_order = atoi(optarg);
                break;
            case 'w':
                window_size = atof(optarg);
                break;
            case 's':
                window_start = atof(optarg);
                break;
            case 'e':
                window_end = atof(optarg);
                break;
            case 'h':
                return EXIT_SUCCESS; 
            default: 
                return EXIT_FAILURE;
        }
    }

    if (file_name.empty()) { 
        cerr << "[ EVALUATOR ]: Please specify a file to process." << endl;
        cerr << "[ EVALUATOR ]: Use --help for usage." << endl;
        return EXIT_FAILURE;
    }

    TFile* file = new TFile(file_name.c_str());
    FlatTupleMaker* tuple = new FlatTupleMaker("fit_results.root", "results"); 

    RootFileReader* reader = new RootFileReader(); 
    reader->parseFile(file);

    BumpHunter* bump_hunter = new BumpHunter(poly_order);
    bump_hunter->setWindowSize(window_size);
    if (bkg_only) bump_hunter->fitBkgOnly();  


    int hist_counter = 0; 
    tuple->addVariable("ap_mass"); 
    tuple->addVariable("yield"); 
    tuple->addVariable("yield_error"); 
    tuple->addVariable("pull");  
    for (auto& hist : reader->get1DHistograms()) { 
        if (hist_counter == hist_count) break;

        map<double, RooFitResult*> results = bump_hunter->fit(hist, window_start, window_end, 0.001);
       
        for (auto& result : results) { 
           
            cout << "Processing result for range: " << result.first << endl;

            RooFitResult* fit_result = result.second; 
            double signal_yield = ((RooRealVar*) fit_result->floatParsFinal().find("bkg yield"))->getVal();
            double signal_yield_error = ((RooRealVar*) fit_result->floatParsFinal().find("bkg yield"))->getError();
           
            tuple->setVariableValue("ap_mass", result.first);  
            tuple->setVariableValue("yield", signal_yield);  
            tuple->setVariableValue("yield_error", signal_yield_error);
            tuple->setVariableValue("pull", signal_yield/signal_yield_error);
              
            tuple->fill(); 

            delete fit_result;  
        }

        results.clear(); 
        hist_counter++;
    }

    tuple->close(); 

    delete bump_hunter; 
    delete reader; 
    delete file;
}
