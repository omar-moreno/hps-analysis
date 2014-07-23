/**
 * @author:	Omar Moreno <omoreno1@ucsc.edu>
 *			Santa Cruz Institute for Particle Physics
 *			University of California, Santa Cruz
 * @date: July 23, 2014
 *
 */

//--- C++ ---//
//-----------//
#include <string>
#include <iostream>

//--- ROOT ---//
//------------//
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>

//--- Utils ---//
//-------------//
#include <PlotUtils.h>

void plotResidualsAtScoringPlane(TFile* root_file, TIter next_key)
{
	TKey* key = NULL; 
	std::string name, class_name;
	TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
	canvas->Print("residuals_at_scoring_plane.pdf[");
	TH1F* histogram = NULL;	
	TF1* fit = NULL;
	while((key=(TKey*) next_key())){

		class_name = key->GetClassName();
		if(class_name.find("TH2") != std::string::npos) continue;

		name = key->GetName();
		if(name.find("scoring") == std::string::npos || name.find("z") != std::string::npos
				|| name.find("_") != std::string::npos){
			continue;
		}

		histogram = (TH1F*) root_file->Get(name.c_str());
		if(name.find("Non") == std::string::npos){
			PlotUtils::set2DPlotStyle(histogram, "x_{extrapolated track position} - x_{MC particle position at scoring plane} (mm)", "Entries");
		} else {
			PlotUtils::set2DPlotStyle(histogram, "y_{extrapolated track position} - y_{MC particle position at scoring plane} (mm)", "Entries");
		}
		histogram->Fit("gaus", "Q", "", histogram->GetMean()-3*histogram->GetRMS(), histogram->GetMean()+3*histogram->GetRMS());
		fit = histogram->GetFunction("gaus");
		fit->SetLineColor(kGreen+3);
		fit->SetLineWidth(2);
		histogram->Draw("E1");
		fit->Draw("same");
		canvas->Print("residuals_at_scoring_plane.pdf(");
	}
	canvas->Print("residuals_at_scoring_plane.pdf]");

	delete canvas;
}

void plotResidualsAtTarget(TFile* root_file, TIter next_key)
{
	TKey* key = NULL; 
	std::string name, class_name;	
	TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
	canvas->Print("residuals_at_target.pdf[");
	TH1F* histogram = NULL;	
	TF1* fit = NULL;
	while((key=(TKey*) next_key())){

		class_name = key->GetClassName(); 
		if(class_name.find("TH2") != std::string::npos) continue;

		name = key->GetName();
		if(name.find("target") == std::string::npos || name.find("z") != std::string::npos
				|| name.find("_") != std::string::npos){
			continue;
		}
		histogram = (TH1F*) root_file->Get(name.c_str());
		if(name.find("Non") == std::string::npos){
			PlotUtils::set2DPlotStyle(histogram, "x_{extrapolated track position} - x_{MC particle position at target} (mm)", "Entries");
		} else {
			PlotUtils::set2DPlotStyle(histogram, "y_{extrapolated track position} - y_{MC particle position at target} (mm)", "Entries");
		}
		histogram->Fit("gaus", "Q", "", histogram->GetMean()-3*histogram->GetRMS(), histogram->GetMean()+3*histogram->GetRMS());
		fit = histogram->GetFunction("gaus");
		fit->SetLineColor(kGreen+3);
		fit->SetLineWidth(2);
		histogram->Draw("E1");
		fit->Draw("same");
		canvas->Print("residuals_at_target.pdf(");
	}
	canvas->Print("residuals_at_target.pdf]");
}

void plotResolutionVsMomentum(TFile* root_file, TIter next_key)
{
	TKey* key = NULL; 
	std::string name, class_name;	
	TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
	canvas->Print("residuals_vs_momentum.pdf[");
	TH1F* histogram = NULL;	
	TF1* fit = NULL;
	TGraphErrors* graph = NULL;
	bool first = true;
	int index = 0;
	std::string graph_name; 
	while((key=(TKey*) next_key())){

		class_name = key->GetClassName(); 
		if(class_name.find("TH2") != std::string::npos) continue;

		name = key->GetName();
		if(name.find("_") == std::string::npos){
			continue;
		}

		index = atoi(name.substr(name.find("_")+1).c_str());

		if(index == 0){
			if(!first){
				graph->Draw("Ap");
				if(name.find("Non") == std::string::npos){
					PlotUtils::set2DPlotStyle(graph, "Track Momentum (GeV)", "#sigma(#Delta x) (mm)");
				} else { 
					PlotUtils::set2DPlotStyle(graph, "Track Momentum (GeV)", "#sigma(#Delta y) (mm)");
				}
				canvas->Print("residuals_vs_momentum.pdf(");
			}
			graph = new TGraphErrors();
			graph_name = name.substr(0, name.find("_"));
			graph->SetTitle(graph_name.c_str()); 
			graph->SetMarkerStyle(20); 
			graph->SetMarkerColor(kBlue-4);
			graph->SetLineWidth(2);
			graph->SetLineColor(kBlue-4);
			first = false; 
		}

		histogram = (TH1F*) root_file->Get(name.c_str());
		histogram->SetTitle((name.substr(0, name.find("_")) + " - projection " + PlotUtils::convertToString(index)).c_str());
		if(name.find("Non") == std::string::npos){
			if(name.find("scoring") != std::string::npos){
				PlotUtils::set2DPlotStyle(histogram, "x_{extrapolated track position} - x_{MC particle position at scoring plane} (mm)", "Entries");
			} else {
				PlotUtils::set2DPlotStyle(histogram, "x_{extrapolated track position} - x_{MC particle position at target} (mm)", "Entries");
			}
		} else {
			if(name.find("scoring") != std::string::npos){
				PlotUtils::set2DPlotStyle(histogram, "y_{extrapolated track position} - y_{MC particle position at scoring plane} (mm)", "Entries");
			} else {
				PlotUtils::set2DPlotStyle(histogram, "y_{extrapolated track position} - y_{MC particle position at target} (mm)", "Entries");
			}
		}
		if(histogram->GetEntries() < 100){
			histogram->Fit("gaus+gaus", "LQ", "", histogram->GetMean()-3*histogram->GetRMS(), histogram->GetMean()+3*histogram->GetRMS());
		} else { 
			histogram->Fit("gaus+gaus", "Q", "", histogram->GetMean()-3*histogram->GetRMS(), histogram->GetMean()+3*histogram->GetRMS());
		}
		fit = histogram->GetFunction("gaus");
		if(histogram->GetEntries() > 10){
			graph->SetPoint(index, .25*index + .25, fit->GetParameter(2));
			graph->SetPointError(index, .25/2, fit->GetParError(2));
		}
		fit->SetLineColor(kGreen+3);
		fit->SetLineWidth(2);
		histogram->Draw("E1");
		fit->Draw("same");
		canvas->Print("residuals_vs_momentum.pdf(");
	}
	graph->Draw("Ap");
	if(name.find("Non") == std::string::npos){
		PlotUtils::set2DPlotStyle(graph, "Track Momentum (GeV)", "#sigma(#Delta x) (mm)");
	} else { 
		PlotUtils::set2DPlotStyle(graph, "Track Momentum (GeV)", "#sigma(#Delta y) (mm)");
	}
	canvas->Print("residuals_vs_momentum.pdf)");
}

void plotAll(std::string root_file_name)
{
	TFile* root_file = new TFile(root_file_name.c_str());
	TIter next_key(root_file->GetListOfKeys());

	PlotUtils::setupStyle();

	plotResidualsAtScoringPlane(root_file, next_key);

	plotResidualsAtTarget(root_file, next_key);

	plotResolutionVsMomentum(root_file, next_key);
}






