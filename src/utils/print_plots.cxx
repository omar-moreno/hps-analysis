/**
 *  @file print_plots.cxx
 *  @brief 
 *  @author Omar Moreno <omoreno1@ucsc.edu>
 *  @date June 27, 2015
 */

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <getopt.h>
#include <iostream>
#include <cstdlib>
#include <list>
#include <fstream>

//--- ROOT ---//
//------------//
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <RootFileReader.h>
#include <PlottingUtils.h>

using namespace std;

int main(int argc, char **argv) {

    string root_file_name;
    string histogram_name; 
    bool save_to_root = false;

    // Parse all the command line arguments. If there are no valid command line
    // arguments passed, print the usage and exit.
    static struct option long_options[] = {
        {"root_file_name", required_argument, 0, 'i' },
        {"name", required_argument, 0, 'n' },
        {"root_file", no_argument, 0, 'r' },
        //{"style", required_argument, 0, 's'},
        //{"pdf", no_argument, 0, 'p' },
        { 0, 0, 0, 0 }
    };

    int option_index = 0;
    int option_char; 
    while ((option_char = getopt_long(argc, argv, "i:n:r", long_options, &option_index)) != -1) {
        switch (option_char) {
            case 'i': 
                root_file_name = optarg;
                break;
            case 'n':
                histogram_name = optarg; 
                break;
            case 'r':
                save_to_root = true;
                break;
            default: 
                break;
        }
    }

    if (root_file_name.empty()) { 
        cerr << "Please specify a ROOT file to process." << endl;
        return EXIT_FAILURE;
    }

    if (histogram_name.empty()) { 
        cerr << "Please specify the name of the histogram you want to print." << endl;
        return EXIT_FAILURE;
    }

    PlottingUtils::setPalette();
    PlottingUtils::setStyle();

    TFile* file = new TFile(root_file_name.c_str());

    RootFileReader* reader = new RootFileReader();
    reader->parseFile(file, histogram_name);

    std::vector<TH1*> histograms_1D; 
    std::vector<TH1*> histograms_2D; 
    std::vector<TGraph*> graphs; 

    histograms_1D = reader->getMatching1DHistograms(histogram_name);
    histograms_2D = reader->getMatching2DHistograms(histogram_name);
    graphs = reader->getMatchingGraphs(histogram_name);

    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
    canvas->Print("plots.pdf[");

    TFile* output_file = NULL;
    if (save_to_root) { 
        output_file = new TFile("root_file.root", "RECREATE");
    }

    for (int hist_n = 0; hist_n < histograms_1D.size(); ++hist_n) { 
        histograms_1D[hist_n]->Draw();
        PlottingUtils::adjust1DPlotRange(histograms_1D[hist_n], 1);
        canvas->Print("plots.pdf(");

        if (save_to_root) histograms_1D[hist_n]->Write();
    }

    for (int hist_n = 0; hist_n < histograms_2D.size(); ++hist_n) { 
        PlottingUtils::adjust2DPlotRange(histograms_2D[hist_n], 1);
        histograms_2D[hist_n]->Draw("colz");
        histograms_2D[hist_n]->GetYaxis()->SetTitleOffset(1.7);
        canvas->Print("plots.pdf(");
        canvas->SaveAs("plot.png");
        if (save_to_root) histograms_2D[hist_n]->Write();
    }

    for (int graph_n = 0; graph_n < graphs.size(); ++graph_n) { 
        graphs[graph_n]->SetName("FEB: 9 Hybrid: 0");
        graphs[graph_n]->SetTitle("FEB: 9 Hybrid: 0");
        graphs[graph_n]->Draw("A*");
        canvas->Print("plots.pdf(");
        if (save_to_root) graphs[graph_n]->Write();
    }



    canvas->Print("plots.pdf]");


    if (save_to_root) output_file->Close();

    return EXIT_SUCCESS;
}
