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

    // Name of file to process
    string file_name;

    // Output file name
    string output_file = "";

    // Default number of histograms to process
    int hist_count = 10;

    // Default polynomial order to use to model the background
    int poly_order = 3;

    // Default mass window size to use when fitting
    double window_size = 0.020; 

    // Default start position of the mass window
    double window_start = 0.03;

    // Maximum position of the end of the mass window
    double window_end = 0.050;

    // Step size
    double step_size = 0.001;

    // Only fit the background
    bool bkg_only = false; 

    // Log fit results
    bool log_fit = false; 

    // Parse all the command line arguments.  If there are no valid command
    // line arguments passed, print the usage and exit the application
    static struct option long_options[] = {
        {"bkg_only",   required_argument, 0, 'b'},
        {"end",        required_argument, 0, 'e'},
        {"file_name",  required_argument, 0, 'i'},
        {"help",       no_argument,       0, 'h'},
        {"log",        no_argument,       0, 'l'},
        {"number",     required_argument, 0, 'n'},
        {"output",     required_argument, 0, 'o'},
        {"poly",       required_argument, 0, 'p'},
        {"start",      required_argument, 0, 's'},
        {"step",       required_argument, 0, 'd'},
        {"window",     required_argument, 0, 'w'},
        {0, 0, 0, 0}
    };
    
    int option_index = 0;
    int option_char; 
    while ((option_char = getopt_long(argc, argv, "be:i:lhn:o:p:s:d:w:", long_options, &option_index)) != -1) {
        switch(option_char) {
            case 'b': 
                bkg_only = true;
                break; 
            case 'e':
                window_end = atof(optarg);
                break;
            case 'i': 
                file_name = optarg;
                break;
            case 'h':
                return EXIT_SUCCESS; 
            case 'l': 
                log_fit = true;
                break; 
            case 'n':
                hist_count = atoi(optarg);
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'p': 
                poly_order = atoi(optarg);
                break;
            case 's':
                window_start = atof(optarg);
                break;
            case 'd':
                step_size = atof(optarg);
                break;
            case 'w':
                window_size = atof(optarg);
                break;
            default: 
                return EXIT_FAILURE;
        }
    }

    // Make sure a file was specified by the user.  If not, warn the user and 
    // exit the application.
    if (file_name.empty()) { 
        cerr << "[ EVALUATOR ]: Please specify a file to process." << endl;
        cerr << "[ EVALUATOR ]: Use --help for usage." << endl;
        return EXIT_FAILURE;
    }

    // Open the ROOT file.  If any problems are encountered, warn the user
    // and exit the application.
    TFile* file = new TFile(file_name.c_str());
    if (file->IsZombie()) { 
        cerr << "[ EVALUATOR ]: Failed to open file " << file_name.c_str(); 
        return EXIT_FAILURE;
    }

    // Create a ROOT file reader and parse the file
    RootFileReader* reader = new RootFileReader(); 
    reader->parseFile(file);

    // Create a new Bump Hunter instance and set the given properties.
    BumpHunter* bump_hunter = new BumpHunter(poly_order);
    bump_hunter->setWindowSize(window_size);
    if (bkg_only) bump_hunter->fitBkgOnly(); 
    if (log_fit) bump_hunter->writeResults();  

    // Build the string that will be used for the results file name
    if (output_file.empty()) { 
        output_file = "order" + to_string(poly_order) + "_window" + to_string(int(window_size*1000)) + "mev"; 
        if (bkg_only) output_file += "_bkg_fit_result.root"; 
        else output_file += "_sig_fit_result.root"; 
    }

    // Create a new flat ntuple and define the variables it will encapsulate.
    FlatTupleMaker* tuple = new FlatTupleMaker(output_file, "results"); 

    tuple->addVariable("hist_n");
    tuple->addVariable("ap_mass"); 
    tuple->addVariable("sig_yield"); 
    tuple->addVariable("sig_yield_error_hi");
    tuple->addVariable("signal_yield_error_low"); 
    tuple->addVariable("bkg_yield"); 
    tuple->addVariable("bkg_yield_error"); 
    tuple->addVariable("nll");
    tuple->addVariable("invalid_nll"); 
    tuple->addVariable("minuit_status");
    tuple->addVariable("edm");
    tuple->addVariable("q0"); 
    tuple->addVariable("p_value"); 
    tuple->addVariable("upper_limit");

    int hist_counter = 0; 
    for (auto& hist : reader->get1DHistograms()) {
        
        // Process only the maximum number of histograms present. 
        if (hist_counter == hist_count) break;

        tuple->setVariableValue("hist_n", (hist_counter + 1)); 
        map<double, HpsFitResult*> results = bump_hunter->fitWindow(hist, window_start, window_end, step_size);
       
        for (auto& result : results) { 
           
            HpsFitResult* fit_result = result.second; 
          
            // Retrieve all of the result of interest. 
            double bkg_yield = ((RooRealVar*) fit_result->getRooFitResult()->floatParsFinal().find("bkg yield"))->getVal();
            double bkg_yield_error = ((RooRealVar*) fit_result->getRooFitResult()->floatParsFinal().find("bkg yield"))->getError();
            double nll = fit_result->getRooFitResult()->minNll();
            double invalid_nll = fit_result->getRooFitResult()->numInvalidNLL();
            double minuit_status = fit_result->getRooFitResult()->status();
            double edm = fit_result->getRooFitResult()->edm(); 

            // Set the values of the results that will be written to the ntuple.
            tuple->setVariableValue("ap_mass", result.first);  
            tuple->setVariableValue("bkg_yield", bkg_yield);  
            tuple->setVariableValue("bkg_yield_error", bkg_yield_error);
            tuple->setVariableValue("nll", nll); 
            tuple->setVariableValue("invalid_nll", invalid_nll); 
            tuple->setVariableValue("minuit_status", minuit_status);
            tuple->setVariableValue("edm", edm);  

            // If this isn't a background only fit evaluation, skip it.
            if (!bkg_only) {
                double signal_yield = ((RooRealVar*) fit_result->getRooFitResult()->floatParsFinal().find("signal yield"))->getVal();
                double signal_yield_error_hi = ((RooRealVar*) fit_result->getRooFitResult()->floatParsFinal().find("signal yield"))->getAsymErrorHi();
                double signal_yield_error_low = ((RooRealVar*) fit_result->getRooFitResult()->floatParsFinal().find("signal yield"))->getAsymErrorLo();
                tuple->setVariableValue("sig_yield", signal_yield);  
                tuple->setVariableValue("sig_yield_error_hi", signal_yield_error_hi);
                tuple->setVariableValue("sig_yield_error_low", signal_yield_error_low);
                tuple->setVariableValue("p_value", fit_result->getPValue());
                tuple->setVariableValue("q0", fit_result->getQ0()); 
                tuple->setVariableValue("upper_limit", fit_result->getUpperLimit());
            }
              
            // Fill the ntuple
            tuple->fill(); 

            // Delete the fit results from memory
            delete fit_result;  
        }

        // Clear the map of results
        results.clear(); 

        // Increment the histogram counter 
        hist_counter++;
    }

    // Close the ntuple
    tuple->close(); 

    delete bump_hunter; 
    delete reader; 
    delete file;
}
