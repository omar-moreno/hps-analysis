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


void ComparePlots::overlayPlots() { 

    TFile* overlay_root_file = new TFile("plot_comparison.root", "RECREATE");
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800); 
    TLegend* legend = NULL;
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
        int max_xbin_value = histogram1D_it->second[0]->GetBinContent(histogram1D_it->second[0]->GetMaximumBin()); 

        for (int hist_n = 1; hist_n < histogram1D_it->second.size(); hist_n++) { 

            histogram1D_it->second[hist_n]->Draw("same");  
            if (histogram1D_it->second[hist_n]->GetBinContent(histogram1D_it->second[hist_n]->GetMaximumBin()) 
                    > max_xbin_value) { 
                max_xbin_value 
                    = histogram1D_it->second[hist_n]->GetBinContent(histogram1D_it->second[hist_n]->GetMaximumBin());
            } 
        }
        histogram1D_it->second[0]->GetYaxis()->SetRangeUser(0, max_xbin_value + .1*max_xbin_value);
        histogram1D_it->second[0]->GetXaxis()->SetRangeUser(2000, 4000);

        if (style.compare("mc") == 0) { 
            
            legend = new TLegend(0.79, 0.79, 0.89, 0.89);
            legend->SetFillColor(0);
            legend->SetLineColor(0);
            legend->AddEntry(histogram1D_it->second[0], "Data", "p");
            legend->AddEntry(histogram1D_it->second[1], "MC", "l");
            legend->Draw();
        
        }

        canvas->Write();
        canvas->Print("plot_comparison.pdf(");
        delete legend; 
    }


    std::map<std::string, std::vector<TH1*> >::iterator histogram2D_it = histogram2D_map.begin();
    for (histogram2D_it; histogram2D_it != histogram2D_map.end(); histogram2D_it++) { 

        histogram2D_it->second[0]->Draw("box");
        int max_bin_value = histogram2D_it->second[0]->GetBinContent(histogram2D_it->second[0]->GetMaximumBin()); 

        for (int hist_n = 1; hist_n < histogram2D_it->second.size(); hist_n++) { 

            histogram2D_it->second[hist_n]->Draw("box same");  
        }
        canvas->Write();
        canvas->Print("plot_comparison.pdf(");
    }


    std::map<std::string, std::vector<TGraph*> >::iterator graph_it = graph_map.begin();
    for (graph_it; graph_it != graph_map.end(); graph_it++) { 

        TMultiGraph* m_graph = new TMultiGraph();
        m_graph->SetNameTitle(graph_it->second[0]->GetName(), graph_it->second[0]->GetTitle());

        TMultiGraph* m_comp = new TMultiGraph();
        m_comp->SetNameTitle(graph_it->second[0]->GetName(), graph_it->second[0]->GetTitle());

        int color_index = 1;
        double x, y, base_x, base_y;
        for (int hist_n = 0; hist_n < graph_it->second.size(); hist_n++) { 

            graph_it->second[hist_n]->SetMarkerStyle(21);
            graph_it->second[hist_n]->SetMarkerColor(color_index);
            graph_it->second[hist_n]->SetMarkerSize(.4);

            m_graph->Add(graph_it->second[hist_n]); 

            TGraph* comparison_graph = new TGraph();
            comparison_graph->SetMarkerStyle(21);
            comparison_graph->SetMarkerColor(color_index);
            comparison_graph->SetLineColor(color_index);
            comparison_graph->SetMarkerSize(.4); 
            color_index++;
            for (int graph_p = 0; graph_p < graph_it->second[hist_n]->GetN(); ++graph_p) { 
               graph_it->second[0]->GetPoint(graph_p, base_x, base_y);
               graph_it->second[hist_n]->GetPoint(graph_p, x, y);
               comparison_graph->SetPoint(graph_p, x, base_y - y);
            }
            m_comp->Add(comparison_graph);
        }
        m_graph->Draw("Ap");
        m_graph->GetXaxis()->SetTitle(graph_it->second[0]->GetXaxis()->GetTitle()); 
        canvas->Write();
        canvas->Print("plot_comparison.pdf(");

        m_comp->Draw("Ap");
        canvas->Print("plot_comparison.pdf(");

        delete m_graph; 
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

    std::map<std::string, std::vector<TGraph*> >::iterator graph_it = graph_map.begin();
    for (graph_it; graph_it != graph_map.end(); graph_it++) { 

        int color_index = 1;
        for (int hist_n = 0; hist_n < graph_it->second.size(); hist_n++) { 

            graph_it->second[hist_n]->SetMarkerStyle(21);
            graph_it->second[hist_n]->SetMarkerColor(color_index);
            graph_it->second[hist_n]->SetLineColor(color_index);
            graph_it->second[hist_n]->SetMarkerSize(.4); 
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
            histogram1D_it->second[hist_n]->SetStats(0);

            histogram1D_it->second[hist_n]->Scale(histogram1D_it->second[0]->Integral()/histogram1D_it->second[hist_n]->Integral());

        }
    }
}
