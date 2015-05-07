/**
 *  @file ComparePlots.h
 *  @brief 
 *  @author Omar Moreno <omoreno1@ucsc.edu>
 *  @date April 29, 2015
 *
 */

#ifndef __COMPARE_PLOTS_H__
#define __COMPARE_PLOTS_H__

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <iostream>
#include <map>
#include <list>
#include <vector>

//------------//
//--- ROOT ---//
//------------//
#include <TCanvas.h>
#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TH1.h>

class ComparePlots { 

    public: 

        /**
         *  Default Constructor
         */
        ComparePlots(); 

        /**
         *  Destructor
         */
        ~ComparePlots();

        /**
         *
         */
        void parseFiles(std::list<TFile*> root_files);

        /**
         *
         */
        void overlayPlots(); 

        /**
         *
         */
        void saveToPdf(std::string pdf_name);

    private: 

        /**
         *
         */
        void addPlots(TList* keys);


        //
        std::map <std::string, std::vector<TH1*> > plot_map;

};

#endif // __COMPARE_PLOTS_H__

