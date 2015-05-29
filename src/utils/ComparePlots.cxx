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
        if (std::string(key->ReadObj()->ClassName()).find("1") == std::string::npos) continue;
        plot_map[key->GetName()].push_back((TH1*) key->ReadObj()); 
        std::cout << "[ ComparePlots ]: Adding file: " << key->GetName() << std::endl;
    } 
}

void ComparePlots::overlayPlots() { 
    
    TFile* overlay_root_file = new TFile("plot_comparison.root", "RECREATE");
    TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500); 
    canvas->Print("plot_comparison.pdf[");
    
    std::map<std::string, std::vector<TH1*> >::iterator plot_it = plot_map.begin();
    for (plot_it; plot_it != plot_map.end(); plot_it++) { 
        
        int color_index = 1;
        std::string options = "";
        if (style.compare("basic") == 0) { 
            plot_it->second[0]->SetLineColor(color_index);
            plot_it->second[0]->SetLineWidth(2);
        } else if (style.compare("mc") == 0) {
            options = "pe"; 
	        plot_it->second[0]->SetMarkerStyle(20);
            plot_it->second[0]->SetMarkerColor(kOrange + 9);
            plot_it->second[0]->SetLineWidth(2);
            plot_it->second[0]->SetLineColor(kOrange + 9);
            plot_it->second[0]->SetMarkerSize(.7); 
        }
        plot_it->second[0]->Draw(options.c_str());
        int max_bin_value = plot_it->second[0]->GetBinContent(plot_it->second[0]->GetMaximumBin()); 

        for (int hist_n = 1; hist_n < plot_it->second.size(); hist_n++) { 
            ++color_index;
            if (style.compare("basic") == 0) { 
                plot_it->second[hist_n]->SetLineColor(color_index);
                plot_it->second[hist_n]->SetLineWidth(2);
            } else if (style.compare("mc") == 0) {
                plot_it->second[hist_n]->SetFillStyle(3003);
                plot_it->second[hist_n]->SetFillColor(kAzure + 2);
                plot_it->second[hist_n]->SetLineColor(kAzure + 3);
                plot_it->second[hist_n]->SetLineWidth(2); 
                plot_it->second[hist_n]->Scale(plot_it->second[0]->Integral()/plot_it->second[hist_n]->Integral());
            }

            plot_it->second[hist_n]->Draw("same");  
            if (plot_it->second[hist_n]->GetBinContent(plot_it->second[hist_n]->GetMaximumBin()) > max_bin_value) { 
                max_bin_value = plot_it->second[hist_n]->GetBinContent(plot_it->second[hist_n]->GetMaximumBin());
            } 
        }
        plot_it->second[0]->GetYaxis()->SetRangeUser(0, max_bin_value + .1*max_bin_value);
        canvas->Write();
        canvas->Print("plot_comparison.pdf(");
    }

    canvas->Print("plot_comparison.pdf]");
    delete overlay_root_file;
    delete canvas; 
}
