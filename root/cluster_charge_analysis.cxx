
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
#include <TGraphErrors.h>

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

RooPlot* processLandauPlot(TH1* landau, double &lm_val, double &le_var);
RooPlot* processSNPlot(TH1* landau, double &lm_val, double &le_var);
std::string getSensorName(std::string histogram_name);

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

    TFile* output_file = new TFile("landau.root", "RECREATE");
    RootFileReader* reader = new RootFileReader(); 
    reader->parseFile(file);

    TCanvas* canvas = new TCanvas("canvas", "canvas", 700, 700);
    TGraphErrors* top_cluster_charge_mpv = new TGraphErrors(18);
    TGraphErrors* bot_cluster_charge_mpv = new TGraphErrors(18);
    top_cluster_charge_mpv->SetTitle("Top SVT Volume");


    bot_cluster_charge_mpv->SetTitle("Bottom SVT Volume");

    //std::vector<TH1*> cluster_amplitude_plots = reader->getMatching1DHistograms("Tracker Cluster Charge");
    //std::vector<TH1*> cluster_amplitude_plots = reader->getMatching1DHistograms("Single Hit Tracker Cluster Charge");
    //std::vector<TH1*> cluster_amplitude_plots = reader->getMatching1DHistograms("Two Hit Tracker Cluster Charge");
    //std::cout << "Size of cluster plots: " << cluster_amplitude_plots.size() << std::endl;

    //std::vector<TH1*>::iterator plot_it = cluster_amplitude_plots.begin();
    canvas->Print("landau_plots.pdf[");
    double lm_val;
    double le_val;
    int top_index = 0;
    int bot_index = 0;
    std::map<int, std::string> top_name_map; 
    std::map<int, std::string> bot_name_map; 
    /*for (plot_it; plot_it != cluster_amplitude_plots.end(); ++plot_it) { 
        
        RooPlot* landau = processLandauPlot(*plot_it, lm_val, le_val);
        std::cout << "Landau mean: " << lm_val << std::endl;
        if (getSensorName((*plot_it)->GetTitle()).find_first_of("t") == 2) { 
            top_cluster_charge_mpv->SetPoint(top_index, top_index, lm_val);
            top_cluster_charge_mpv->SetPointError(top_index, 0, le_val);
            top_name_map[top_index] = getSensorName((*plot_it)->GetTitle());
            top_index++; 
        } else { 
            bot_cluster_charge_mpv->SetPoint(bot_index, bot_index, lm_val);
            bot_cluster_charge_mpv->SetPointError(bot_index, 0, le_val);
            bot_name_map[bot_index] = getSensorName((*plot_it)->GetTitle());
            bot_index++; 
        }
        landau->Draw("");
        canvas->Print("landau_plots.pdf(");
        landau->Write();
        delete landau;
    }

    for (int index = 0; index < 18; index++) { 
        if (index == 0) continue; 
        int label_bin = top_cluster_charge_mpv->GetXaxis()->FindBin(index);
        top_cluster_charge_mpv->GetXaxis()->SetBinLabel(label_bin, top_name_map[index].c_str() );
        top_cluster_charge_mpv->GetXaxis()->SetLabelSize(0.02);

        label_bin = bot_cluster_charge_mpv->GetXaxis()->FindBin(index);
        bot_cluster_charge_mpv->GetXaxis()->SetBinLabel(label_bin, bot_name_map[index].c_str() );
        bot_cluster_charge_mpv->GetXaxis()->SetLabelSize(0.02);
    }


    top_cluster_charge_mpv->Draw("A*");
    top_cluster_charge_mpv->GetYaxis()->SetTitle("Cluster Amplitude (ADC Counts)");
    top_cluster_charge_mpv->Write();
    canvas->Print("landau_plots.pdf(");


    bot_cluster_charge_mpv->Draw("A*");
    bot_cluster_charge_mpv->GetYaxis()->SetTitle("Cluster Amplitude (ADC Counts)");
    bot_cluster_charge_mpv->Write();
    canvas->Print("landau_plots.pdf(");


    delete top_cluster_charge_mpv;
*/

    std::vector<TH1*> sn_plots = reader->getMatching1DHistograms("Track Two Hit Signal to Noise");
            
    //std::vector<TH1*>::iterator plot_it = cluster_amplitude_plots.begin();
    std::cout << "Size of sn plots: " << sn_plots.size() << std::endl;
    std::vector<TH1*>::iterator plot_it = sn_plots.begin();

    TGraphErrors* top_cluster_signal_mpv = new TGraphErrors(18);
    TGraphErrors* bot_cluster_signal_mpv = new TGraphErrors(18);
    top_cluster_charge_mpv->SetTitle("Top SVT Volume");


    bot_cluster_charge_mpv->SetTitle("Bottom SVT Volume");


    for (plot_it; plot_it != sn_plots.end(); ++plot_it) { 
        
        RooPlot* landau = processSNPlot(*plot_it, lm_val, le_val);
        std::cout << "Landau mean: " << lm_val << std::endl;
        if (getSensorName((*plot_it)->GetTitle()).find_first_of("t") == 2) { 
            top_cluster_signal_mpv->SetPoint(top_index, top_index, lm_val);
            top_cluster_signal_mpv->SetPointError(top_index, 0, le_val);
            top_name_map[top_index] = getSensorName((*plot_it)->GetTitle());
            top_index++; 
        } else { 
            bot_cluster_signal_mpv->SetPoint(bot_index, bot_index, lm_val);
            bot_cluster_signal_mpv->SetPointError(bot_index, 0, le_val);
            bot_name_map[bot_index] = getSensorName((*plot_it)->GetTitle());
            bot_index++; 
        }
        landau->Draw("");
        canvas->Print("landau_plots.pdf(");
        landau->Write();
        delete landau;
    }


    for (int index = 0; index < 18; index++) { 
        if (index == 0) continue; 
        int label_bin = top_cluster_signal_mpv->GetXaxis()->FindBin(index);
        top_cluster_signal_mpv->GetXaxis()->SetBinLabel(label_bin, top_name_map[index].c_str() );
        top_cluster_signal_mpv->GetXaxis()->SetLabelSize(0.02);

        label_bin = bot_cluster_signal_mpv->GetXaxis()->FindBin(index);
        bot_cluster_signal_mpv->GetXaxis()->SetBinLabel(label_bin, bot_name_map[index].c_str() );
        bot_cluster_signal_mpv->GetXaxis()->SetLabelSize(0.02);
    }

    top_cluster_signal_mpv->Draw("A*");
    top_cluster_signal_mpv->GetYaxis()->SetTitle("Signal To Noise");
    top_cluster_signal_mpv->Write();
    canvas->Print("landau_plots.pdf(");


    bot_cluster_signal_mpv->Draw("A*");
    bot_cluster_signal_mpv->GetYaxis()->SetTitle("Signal to Noise");
    bot_cluster_signal_mpv->Write();
    canvas->Print("landau_plots.pdf(");


    canvas->Print("landau_plots.pdf]");
    output_file->Close();
    file->Close();
    delete reader;

    return EXIT_SUCCESS;
}


RooPlot* processLandauPlot(TH1* landau, double &lm_val, double &le_val) {

    
    RooRealVar cluster_charge("cluster_charge", "Cluster Amplitude (ADC Counts)", 0, 5000);
    RooDataHist landau_data("Cluster Charge", "Cluster Charge", cluster_charge, landau); 

    RooRealVar landau_mean("landau_mean", "landau_mean", 1500., 1000, 2500);
    RooRealVar landau_sigma("landau_sigma", "landau_sigma", 250, 10, 500);
    RooLandau landau_pdf("landau_pdf", "landau_pdf", cluster_charge, landau_mean, landau_sigma);

    RooRealVar gauss_mean("gauss_mean", "gauss_mean", 0 );
    RooRealVar gauss_sigma("gauss_sigma", "gauss_sigma", 60, 10, 400);
    RooGaussian gauss("gauss", "gauss", cluster_charge, gauss_mean, gauss_sigma);

    RooNumConvPdf landau_x_gauss("landau_x_gauss", "landau_x_gauss", cluster_charge, landau_pdf, gauss);
    landau_x_gauss.setConvolutionWindow(landau_mean, landau_sigma, 17);

    RooFitResult* result = landau_x_gauss.fitTo(landau_data, Range(1200., 3000.), Save());

    RooPlot* cluster_charge_frame = cluster_charge.frame();
    cluster_charge_frame->SetTitle(landau->GetTitle());
    cluster_charge_frame->SetTitle(getSensorName(landau->GetTitle()).c_str());
    //cluster_charge_frame->SetName(landau->GetName());
    landau_data.plotOn(cluster_charge_frame);
    landau_pdf.plotOn(cluster_charge_frame, LineStyle(kDashed), LineColor(kRed));
    //gauss.plotOn(cluster_charge_frame, LineStyle(kDashed));
    landau_x_gauss.plotOn(cluster_charge_frame, Range(0., 5000.));
   
    RooRealVar* mean = (RooRealVar*) result->floatParsFinal().find("landau_mean");
    RooRealVar* err = (RooRealVar*) mean->errorVar();
    std::cout << "Landau mean: " << mean->getVal() << std::endl;
    lm_val = mean->getVal();
    le_val = err->getVal();

    delete err;

    //    fit_params.Print(); 
    return cluster_charge_frame;
}


RooPlot* processSNPlot(TH1* landau, double &lm_val, double &le_val) {

    
    RooRealVar signal_to_noise("signal_to_noise", "Signal to Noise", 0, 50);
    RooDataHist landau_data("Signal To Noise", "Signal To Noise", signal_to_noise, landau); 

    RooRealVar landau_mean("landau_mean", "landau_mean", 25, 15, 30);
    RooRealVar landau_sigma("landau_sigma", "landau_sigma", 5, .1, 15);
    RooLandau landau_pdf("landau_pdf", "landau_pdf", signal_to_noise, landau_mean, landau_sigma);

    RooRealVar gauss_mean("gauss_mean", "gauss_mean", 0 );
    RooRealVar gauss_sigma("gauss_sigma", "gauss_sigma", 1, 0.1, 20);
    RooGaussian gauss("gauss", "gauss", signal_to_noise, gauss_mean, gauss_sigma);

    RooNumConvPdf landau_x_gauss("landau_x_gauss", "landau_x_gauss", signal_to_noise, landau_pdf, gauss);
    landau_x_gauss.setConvolutionWindow(landau_mean, landau_sigma, 12);

    RooFitResult* result = landau_x_gauss.fitTo(landau_data, Range(12, 40), Save());

    RooPlot* signal_to_noise_frame = signal_to_noise.frame();
    signal_to_noise_frame->SetTitle(landau->GetTitle());
    signal_to_noise_frame->SetTitle(getSensorName(landau->GetTitle()).c_str());
    //signal_to_noise_frame->SetName(landau->GetName());
    landau_data.plotOn(signal_to_noise_frame);
    landau_pdf.plotOn(signal_to_noise_frame, LineStyle(kDashed), LineColor(kRed));
    //gauss.plotOn(signal_to_noise_frame, LineStyle(kDashed));
    landau_x_gauss.plotOn(signal_to_noise_frame, Range(0., 50.));
   
    RooRealVar* mean = (RooRealVar*) result->floatParsFinal().find("landau_mean");
    RooRealVar* err = (RooRealVar*) mean->errorVar();
    std::cout << "Landau mean: " << mean->getVal() << std::endl;
    lm_val = mean->getVal();
    le_val = err->getVal();

    delete err;

    //    fit_params.Print(); 
    return signal_to_noise_frame;
}



std::string getSensorName(std::string histogram_name) { 
    
    std::string buffer = histogram_name.substr(7, 3) + " ";
    
    std::size_t find_index;  
    if (histogram_name.find("axial") != std::string::npos) { 
        find_index = histogram_name.find("axial");
        buffer += histogram_name.substr(find_index, 5); 
    } else if (histogram_name.find("stereo") != std::string::npos) { 
        find_index = histogram_name.find("stereo");
        buffer += histogram_name.substr(find_index, 6); 
    }

    buffer += " ";
    if (histogram_name.find("hole") != std::string::npos) { 
        find_index = histogram_name.find("hole");
        buffer += histogram_name.substr(find_index, 4); 
    } else if (histogram_name.find("slot") != std::string::npos) { 
        find_index = histogram_name.find("slot");
        buffer += histogram_name.substr(find_index, 4); 
    }

    return buffer;
}

