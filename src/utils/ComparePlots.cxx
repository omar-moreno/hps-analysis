/**
 *  @file ComparePlots.cxx
 *  @brief 
 *  @author Omar Moreno <omoreno1@ucsc.edu>
 *  @date April 29, 2015
 *
 */

#include <ComparePlots.h>

ComparePlots::ComparePlots() { 
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
        plot_map[key->GetName()].push_back((TH1*) key->ReadObj()); 
    } 
}

void ComparePlots::overlayPlots() { 
    
    TFile* overlay_root_file = new TFile("plot_comparison.root", "RECREATE");
    TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500); 
    
    std::map<std::string, std::vector<TH1*> >::iterator plot_it = plot_map.begin();
    for (plot_it; plot_it != plot_map.end(); plot_it++) { 
        int color_index = 1;
        plot_it->second[0]->SetLineColor(color_index);
	    plot_it->second[0]->Draw();
        for (int hist_n = 1; hist_n < plot_it->second.size(); hist_n++) { 
            ++color_index;
            plot_it->second[hist_n]->SetLineColor(color_index);
            plot_it->second[hist_n]->Draw("same");  
        }
        canvas->Write();
    }

    overlay_root_file->Close();

    delete overlay_root_file;
    delete canvas; 
}

void ComparePlots::saveToPdf(std::string pdf_name) { 
   
    TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500); 
    std::map<std::string, std::vector<TH1*> >::iterator plot_it = plot_map.begin();
    canvas->Print((pdf_name + "[").c_str());
    for (plot_it; plot_it != plot_map.end(); ++plot_it) { 
        for (int hist_n = 0; hist_n < plot_it->second.size(); ++hist_n) { 
            plot_it->second[hist_n]->Draw();
            canvas->Print((pdf_name + "(").c_str());
        }
    }
    canvas->Print((pdf_name + "]").c_str());
    delete canvas;
}
