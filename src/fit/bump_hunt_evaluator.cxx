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
#include <Plotter.h>

using namespace std;

int main(int argc, char **argv) { 

    string file_name;
    int hist_count = 10;
    int poly_order = 3;
    double window_size = 0.020; 
    double window_start = 0.03;
    double window_end = 0.050;

    // Parse all the command line arguments.  If there are no valid command
    // line arguments passed, print the usage and exit the application
    static struct option long_options[] = {
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
    while ((option_char = getopt_long(argc, argv, "i:n:o:w:s:e:h", long_options, &option_index)) != -1) {

        switch(option_char){
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

    RootFileReader* reader = new RootFileReader(); 
    reader->parseFile(file);

    BumpHunter* bump_hunter = new BumpHunter(poly_order);
    bump_hunter->setWindowSize(window_size); 

    Plotter* plotter = new Plotter(); 

    int hist_counter = 0; 
    string hist_name; 
    for (auto& hist : reader->get1DHistograms()) { 
        cout << "Histogram: " << hist->GetName() << endl;
        if (hist_counter == hist_count) break;

        map<string, RooFitResult*> results = bump_hunter->fit(hist, window_start, window_end, 0.01);
       
        for (auto& result : results) { 
            
            hist_name = result.first + " - Signal Yield"; 
            cout << "Processing result for range: " << result.first << endl;
            if (!plotter->has1DHistogram(hist_name)) { 
                plotter->build1DHistogram(hist_name, 100, -3000, 3000); 
                plotter->build1DHistogram(result.first + " - Signal Yield Error", 50, -1000, 1000); 
                plotter->build1DHistogram(result.first + " - Pull", 50, -10, 10); 
            }

            RooFitResult* fit_result = result.second; 
            double signal_yield = ((RooRealVar*) fit_result->floatParsFinal().find("signal yield"))->getVal();
            double signal_yield_error = ((RooRealVar*) fit_result->floatParsFinal().find("signal yield"))->getError();
            plotter->get1DHistogram(result.first + " - Signal Yield")->Fill(signal_yield);
            plotter->get1DHistogram(result.first + " - Signal Yield Error")->Fill(signal_yield_error);
            plotter->get1DHistogram(result.first + " - Pull")->Fill(signal_yield/signal_yield_error);
            delete fit_result;  
        }

        results.clear(); 
        hist_counter++;
    }

    plotter->saveToRootFile("fit_results.root"); 

    delete plotter; 
    delete bump_hunter; 
    delete reader; 
    delete file;
}
