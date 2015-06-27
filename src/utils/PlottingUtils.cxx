/**
 *	@file	PlottingUtils.h
 *	@brief	A compilation of plotting utilities.
 *	@author Omar Moreno <omoreno1@ucsc.edu>
 *			Santa Cruz Institute for Paritcle Physics
 *			University of California, Santa Cruz
 *	@date	2013
 *
 */

#include <PlottingUtils.h>

void PlottingUtils::setPalette() {
        
		// Setup the color palette
        const int n_rgb  = 5; 
        const int n_cont = 255;

        double stops[n_rgb] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
        double red[n_rgb]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
        double green[n_rgb] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
        double blue[n_rgb]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
        TColor::CreateGradientColorTable(n_rgb, stops, red, green, blue, n_cont);
        gStyle->SetNumberContours(n_cont); 
}

void PlottingUtils::setStyle() { 
	
	// Set frame attributes
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameBorderSize(0);
	gStyle->SetFrameFillColor(0);

	// Set canvas attributes
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasBorderSize(0);
	gStyle->SetCanvasColor(0);

	// Set axis attributes
	int font = 42;
	double font_size = 0.03;
	double offset = 1.29;
	gStyle->SetLabelFont(font, "x");
	gStyle->SetLabelSize(font_size, "x"); 
	gStyle->SetTitleFont(font, "x");
	gStyle->SetTitleSize(font_size, "x");
	gStyle->SetTitleOffset(offset, "x");	
	gStyle->SetLabelFont(font, "y");
	gStyle->SetLabelSize(font_size, "y"); 
	gStyle->SetTitleFont(font, "y");
	gStyle->SetTitleSize(font_size, "y"); 
	gStyle->SetTitleOffset(offset + .50, "y");	
	gStyle->SetLabelFont(font, "z");
	gStyle->SetTitleFont(font, "y");

	// Enable tick marks on the top and right of plots
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1); 

	// Set marker attributes
	gStyle->SetMarkerStyle(20); 
	gStyle->SetMarkerSize(.15);

	gStyle->SetOptStat(0);

	// Set margin attributes
	gStyle->SetPadLeftMargin(.12);	
	gStyle->SetPadRightMargin(.12);	
}	
