/**
 *  @file ComparePlots.cxx
 *  @brief 
 *  @author Omar Moreno <omoreno1@ucsc.edu>
 *  @date April 29, 2015
 *
 */

#include <ComparePlots.h>

ComparePlots::ComparePlots() 
    : style("") { 
    }

ComparePlots::~ComparePlots() { 
}

void ComparePlots::parseFiles(std::list<TFile*> root_files) { 

    // Loop over all of the ROOT files and create the plot maps
    std::list<TFile*>::iterator root_files_it = root_files.begin();
    for (root_files_it; root_files_it != root_files.end(); ++root_files_it) { 
        std::cout << "[ ComparePlots ]: Processing file: " << (*root_files_it)->GetName() << std::endl;        
        this->addPlots((*root_files_it)->GetListOfKeys());    
    }
}

void ComparePlots::addPlots(TList* keys) {

    TIter next(keys);
    while (TKey *key = (TKey*) next()) { 
        if (key->IsFolder()) this->addPlots(((TDirectory*) key->ReadObj())->GetListOfKeys());
        if (std::string(key->ReadObj()->ClassName()).find("1") != std::string::npos) {
            histogram1D_map[key->GetName()].push_back((TH1*) key->ReadObj()); 
        } else if (std::string(key->ReadObj()->ClassName()).find("2") != std::string::npos) {
            histogram2D_map[key->GetName()].push_back((TH1*) key->ReadObj()); 
        }
        std::cout << "[ ComparePlots ]: Adding file: " << key->GetName() << std::endl;
    } 
}

void ComparePlots::overlayPlots() { 

    TFile* overlay_root_file = new TFile("plot_comparison.root", "RECREATE");
    TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500); 
    canvas->Print("plot_comparison.pdf[");
    
    std::string options = ""; 
    this->applyBasic2DStyle();
    if (style.compare("basic") == 0) { 
        this->applyBasic1DStyle();
    } else if (style.compare("mc") == 0) {
        options = "pe"; 
        this->applyMCStyle();
    }

    std::map<std::string, std::vector<TH1*> >::iterator histogram1D_it = histogram1D_map.begin();
    for (histogram1D_it; histogram1D_it != histogram1D_map.end(); histogram1D_it++) { 

        histogram1D_it->second[0]->Draw(options.c_str());
        int max_bin_value = histogram1D_it->second[0]->GetBinContent(histogram1D_it->second[0]->GetMaximumBin()); 

        for (int hist_n = 1; hist_n < histogram1D_it->second.size(); hist_n++) { 

            histogram1D_it->second[hist_n]->Draw("same");  
            if (histogram1D_it->second[hist_n]->GetBinContent(histogram1D_it->second[hist_n]->GetMaximumBin()) 
                    > max_bin_value) { 
                max_bin_value 
                    = histogram1D_it->second[hist_n]->GetBinContent(histogram1D_it->second[hist_n]->GetMaximumBin());
            } 
        }
        histogram1D_it->second[0]->GetYaxis()->SetRangeUser(0, max_bin_value + .1*max_bin_value);
        canvas->Write();
        canvas->Print("plot_comparison.pdf(");
    }


    std::map<std::string, std::vector<TH1*> >::iterator histogram2D_it = histogram2D_map.begin();
    for (histogram2D_it; histogram2D_it != histogram2D_map.end(); histogram2D_it++) { 

        histogram2D_it->second[0]->Draw("box");
        int max_bin_value = histogram2D_it->second[0]->GetBinContent(histogram2D_it->second[0]->GetMaximumBin()); 

        for (int hist_n = 1; hist_n < histogram2D_it->second.size(); hist_n++) { 

            histogram2D_it->second[hist_n]->Draw("box same");  
            /*if (histogram2D_it->second[hist_n]->GetBinContent(histogram2D_it->second[hist_n]->GetMaximumBin()) 
                    > max_bin_value) { 
                max_bin_value 
                    = histogram2D_it->second[hist_n]->GetBinContent(histogram2D_it->second[hist_n]->GetMaximumBin());
            }*/ 
        }
        //histogram2D_it->second[0]->GetYaxis()->SetRangeUser(0, max_bin_value + .1*max_bin_value);
        canvas->Write();
        canvas->Print("plot_comparison.pdf(");
    }


    canvas->Print("plot_comparison.pdf]");
    delete overlay_root_file;
    delete canvas; 
}

void ComparePlots::applyBasic1DStyle() {

    std::map<std::string, std::vector<TH1*> >::iterator histogram1D_it = histogram1D_map.begin();
    for (histogram1D_it; histogram1D_it != histogram1D_map.end(); histogram1D_it++) { 

        int color_index = 1;
        for (int hist_n = 0; hist_n < histogram1D_it->second.size(); ++hist_n) {
            histogram1D_it->second[hist_n]->SetLineColor(color_index);
            histogram1D_it->second[hist_n]->SetFillStyle(3003);
            histogram1D_it->second[hist_n]->SetFillColor(color_index);
            histogram1D_it->second[hist_n]->SetLineWidth(2);
            color_index++;
        }
    }
}

void ComparePlots::applyBasic2DStyle() { 

    std::map<std::string, std::vector<TH1*> >::iterator histogram2D_it = histogram2D_map.begin();
    for (histogram2D_it; histogram2D_it != histogram2D_map.end(); histogram2D_it++) { 

        int color_index = 1;
        for (int hist_n = 0; hist_n < histogram2D_it->second.size(); ++hist_n) {
            histogram2D_it->second[hist_n]->SetLineColor(color_index);
            histogram2D_it->second[hist_n]->SetLineWidth(2);
            color_index++;
        }
    }
}

void ComparePlots::applyMCStyle() {

    std::map<std::string, std::vector<TH1*> >::iterator histogram1D_it = histogram1D_map.begin();
    for (histogram1D_it; histogram1D_it != histogram1D_map.end(); histogram1D_it++) { 

        int color_index = 0;
        for (int hist_n = 0; hist_n < histogram1D_it->second.size(); ++hist_n) {
            if (hist_n == 0) {
                color_index = kOrange + 9;
            } else { 
                color_index = kAzure + 3;
            }
            histogram1D_it->second[hist_n]->SetMarkerStyle(20);
            histogram1D_it->second[hist_n]->SetMarkerColor(color_index);
            histogram1D_it->second[hist_n]->SetLineWidth(2);
            histogram1D_it->second[hist_n]->SetLineColor(color_index);
            histogram1D_it->second[hist_n]->SetMarkerSize(.7); 
            histogram1D_it->second[hist_n]->SetFillStyle(3003);
            histogram1D_it->second[hist_n]->SetFillColor(color_index - 1);
            
            histogram1D_it->second[hist_n]->Scale(histogram1D_it->second[0]->Integral()/histogram1D_it->second[hist_n]->Integral());

        }
    }
}
