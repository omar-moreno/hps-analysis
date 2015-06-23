
//------------------//
//--- C++ StdLib ---//
//------------------//
#include <cstdlib>
#include <getopt.h>

//------------//
//--- ROOT ---//
//------------//
#include <TFile.h>
#include <TCanvas.h>

//--------------//
//--- RooFit ---//
//--------------//
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooLandau.h>
#include <RooGaussian.h>
#include <RooNumConvPdf.h>
#include <RooAbsReal.h>

//-------------//
//--- Utils ---//
//-------------//
#include <RootFileReader.h>

using namespace std;
using namespace RooFit;

RooPlot* processLandauPlot(TH1* landau);

int main (int argc, char** argv) { 
  
    string file_name;

    // Parse all the command line arguments.  If there are no valid command
    // line arguments passed, print the usage and exit the application
    static struct option long_options[] = {
        {"file_name",  required_argument, 0, 'i'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int option_char; 
    while ((option_char = getopt_long(argc, argv, "i:", long_options, &option_index)) != -1) {

        switch(option_char) {
            case 'i': 
                file_name = optarg;
                break;
            default: 
                return EXIT_FAILURE;
        }
    }

    // Open the ROOT file.  If the file can't be opened, exit the 
    // application
    TFile* file = new TFile(file_name.c_str());

    RootFileReader* reader = new RootFileReader(); 
    reader->parseFile(file);

    TCanvas* canvas = new TCanvas("canvas", "canvas", 700, 700);

    std::vector<TH1*> cluster_amplitude_plots = reader->getMatching1DHistograms("Tracker Cluster Charge");
    std::cout << "Size of cluster plots: " << cluster_amplitude_plots.size() << std::endl;

    std::vector<TH1*>::iterator plot_it = cluster_amplitude_plots.begin();
    canvas->Print("landau_plots.pdf[");
    for (plot_it; plot_it != cluster_amplitude_plots.end(); ++plot_it) { 
        
        RooPlot* landau = processLandauPlot(*plot_it);
        landau->Draw();
        canvas->Print("landau_plots.pdf(");
        delete landau;
    }
    canvas->Print("landau_plots.pdf]");



    file->Close();
    delete reader;

    return EXIT_SUCCESS;
}


RooPlot* processLandauPlot(TH1* landau) {

    
    RooRealVar cluster_charge("cluster_charge", "cluster_charge", 0, 5000);
    RooDataHist landau_data("Cluster Charge", "Cluster Charge", cluster_charge, landau); 

    RooRealVar landau_mean("landau_mean", "landau_mean", 1500., 1000, 2500);
    RooRealVar landau_sigma("landau_sigma", "landau_sigma", 250, 10, 500);
    RooLandau landau_pdf("landau", "landau", cluster_charge, landau_mean, landau_sigma);
    //landau_pdf.setConvolutionWindow(1500, 500, 5);

    RooRealVar gauss_mean("gauss_mean", "gauss_mean", 0 );
    RooRealVar gauss_sigma("gauss_sigma", "gauss_sigma", 60, 10, 200);
    RooGaussian gauss("gauss", "gauss", cluster_charge, gauss_mean, gauss_sigma);
    //gauss.setConvolutionWindow(1500, 500, 5);

    RooNumConvPdf lxg("lxg", "lxg", cluster_charge, landau_pdf, gauss);
    lxg.setConvolutionWindow(landau_mean, landau_sigma, 20);

    lxg.fitTo(landau_data, Range(1000, 5000));
    //landau_pdf.fitTo(landau_data, Range(1000, 5000));

    RooPlot* cluster_charge_frame = cluster_charge.frame();
    landau_data.plotOn(cluster_charge_frame);
    landau_pdf.plotOn(cluster_charge_frame, LineStyle(kDashed), LineColor(kRed));
    //gauss.plotOn(cluster_charge_frame, LineStyle(kDashed));
    lxg.plotOn(cluster_charge_frame);

    return cluster_charge_frame;
}

