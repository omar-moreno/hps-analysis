
#ifndef __PLOTTER_H__
#define __PLOTTER_H__

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <iostream>
#include <typeinfo>
#include <cstdlib>

//------------//
//--- ROOT ---//
//------------//
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TCanvas.h>

class Plotter { 

    public: 

        /**
         * Constructor
         */
        Plotter();

        /**
         *
         */
        ~Plotter(); 

        /**
         *
         */
        Plotter* setType(std::string type) { this->type = type; return this; };

        /**
         *
         */
        Plotter* setLineColor(int color) { this->color = color; return this; };

        /**
         *
         */
        TH1* build1DHistogram(std::string name, int n_bins, double x_min, double x_max); 
       
        /**
         *
         */ 
        TH2* build2DHistogram(std::string name, int n_bins_x, double x_min, double x_max,
                int n_bins_y, double y_min, double y_max);
        /**
         *
         */
        TH1* get1DHistogram(std::string name);
        
        /**
         *
         */
        TH2* get2DHistogram(std::string name);

        /**
         *
         */
        void saveToRootFile(std::string file_name); 
        
        /**
         *
         */
        void saveToPdf(std::string file_name); 



    private: 

        std::unordered_map<std::string, TH1*> histogram1D_map;
        std::unordered_map<std::string, TH2*> histogram2D_map;

        std::string type; 

        int color; 

}; // Plotter

#endif // __PLOTTER_H__
