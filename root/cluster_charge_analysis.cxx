
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
#include <RooFitResult.h>
#include <RooArgList.h>

//-------------//
//--- Utils ---//
//-------------//
#include <RootFileReader.h>

using namespace std;
using namespace RooFit;

RooPlot* processLandauPlot(TH1* landau, double &lm_val);

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
    TGraph* top_cluster_charge_mpv = new TGraph(36);

    std::vector<TH1*> cluster_amplitude_plots = reader->getMatching1DHistograms("Tracker Cluster Charge");
    std::cout << "Size of cluster plots: " << cluster_amplitude_plots.size() << std::endl;

    std::vector<TH1*>::iterator plot_it = cluster_amplitude_plots.begin();
    canvas->Print("landau_plots.pdf[");
    double lm_val;
    int index = 0; 
    for (plot_it; plot_it != cluster_amplitude_plots.end(); ++plot_it) { 
        
        RooPlot* landau = processLandauPlot(*plot_it, lm_val);
        std::cout << "Landau mean: " << lm_val << std::endl;
        top_cluster_charge_mpv->SetPoint(index, index, lm_val);
        index++;
        landau->Draw("");
        canvas->Print("landau_plots.pdf(");
        delete landau;
        if (index == 5) break;
    }
    top_cluster_charge_mpv->Draw("A*");
    canvas->Print("landau_plots.pdf(");
    canvas->Print("landau_plots.pdf]");

    file->Close();
    delete reader;
    delete top_cluster_charge_mpv;

    return EXIT_SUCCESS;
}


RooPlot* processLandauPlot(TH1* landau, double &lm_val) {

    
    RooRealVar cluster_charge("cluster_charge", "cluster_charge", 0, 5000);
    RooDataHist landau_data("Cluster Charge", "Cluster Charge", cluster_charge, landau); 

    RooRealVar landau_mean("landau_mean", "landau_mean", 1500., 1000, 2500);
    RooRealVar landau_sigma("landau_sigma", "landau_sigma", 250, 10, 500);
    RooLandau landau_pdf("landau_pdf", "landau_pdf", cluster_charge, landau_mean, landau_sigma);

    RooRealVar gauss_mean("gauss_mean", "gauss_mean", 0 );
    RooRealVar gauss_sigma("gauss_sigma", "gauss_sigma", 60, 10, 300);
    RooGaussian gauss("gauss", "gauss", cluster_charge, gauss_mean, gauss_sigma);

    RooNumConvPdf landau_x_gauss("landau_x_gauss", "landau_x_gauss", cluster_charge, landau_pdf, gauss);
    landau_x_gauss.setConvolutionWindow(landau_mean, landau_sigma, 20);

    RooFitResult* result = landau_x_gauss.fitTo(landau_data, Range(1000., 5000.), Save());

    RooPlot* cluster_charge_frame = cluster_charge.frame();
    landau_data.plotOn(cluster_charge_frame);
    landau_pdf.plotOn(cluster_charge_frame, LineStyle(kDashed), LineColor(kRed));
    gauss.plotOn(cluster_charge_frame, LineStyle(kDashed));
    landau_x_gauss.plotOn(cluster_charge_frame, Range(0., 5000.));
   
    RooRealVar* mean = (RooRealVar*) result->floatParsFinal().find("landau_mean");
    std::cout << "Landau mean: " << mean->getVal() << std::endl;
    lm_val = mean->getVal();

    //    fit_params.Print(); 
    return cluster_charge_frame;
}

