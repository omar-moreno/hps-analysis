
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
#include <RooGaussian.h>
#include <RooNumConvPdf.h>
#include <RooAbsReal.h>
#include <RooFitResult.h>
#include <RooArgList.h>
#include <RooAddPdf.h>

using namespace std;
using namespace RooFit;

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

    TCanvas* canvas = new TCanvas("canvas", "canvas", 700, 700);
    canvas->Print("gbl_fits.pdf["); 

    TH1* p_histo = (TH1*) file->Get("p - bottom");

    RooRealVar p("fee_p", "Momentum (GeV)", 0, 2);
    RooDataHist p_data("fee_p_data", "FEE Momentum", p, p_histo); 
    
    RooRealVar sig_gauss_mean("sig_gauss_mean", "sig_gauss_mean", .99, .8, 1.3);
    RooRealVar sig_gauss_sigma("sig_gauss_sigma", "sig_gauss_sigma", .05, 0, 2);
    RooGaussian sig_gauss("sig_gauss", "sig_gauss", p, sig_gauss_mean, sig_gauss_sigma);

    RooRealVar bkg_gauss_mean("bkg_gauss_mean", "bkg_gauss_mean", .5, 0, .7);
    RooRealVar bkg_gauss_sigma("bkg_gauss_sigma", "bkg_gauss_sigma", 5, 0, 100);
    RooGaussian bkg_gauss("bkg_gauss", "bkg_gauss", p, bkg_gauss_mean, bkg_gauss_sigma);
   
    //RooRealVar bgk_argus_mean(
    //RooArgusBG("bkg_argus", "bkg_argus", p, bgk_gauss_mean, bkg_

    RooRealVar n_sig("n_sig", "#signal events", 500., 0, 100000);
    RooRealVar n_bkg("n_bkg", "#background events", 500, 0, 10000);  
    RooAddPdf sum("sum", "g+g", RooArgList(sig_gauss, bkg_gauss), RooArgList(n_sig, n_bkg)); 

    RooNumConvPdf gauss_x_gauss("gauss_x_gauss", "gauss_x_gauss", p, sig_gauss, bkg_gauss);
    //gauss_x_gauss.setConvolutionWindow(sig_gauss_mean, sig_gauss_sigma, 10);

    //gauss_x_gauss.fitTo(p_data,Save());  
    sum.fitTo(p_data,Save());  

    //sig_gauss.fitTo(p_data, Range(.6, 1.5), Save());

    RooPlot* p_frame = p.frame();
    p_data.plotOn(p_frame); 
    //sig_gauss.plotOn(p_frame); 
    //gauss_x_gauss.plotOn(p_frame); 
    sum.plotOn(p_frame);

    p_frame->Draw(""); 
    
    canvas->Print("gbl_fits.pdf("); 
    canvas->Print("gbl_fits.pdf]");

    file->Close();
    
   return EXIT_SUCCESS;  
    
}
