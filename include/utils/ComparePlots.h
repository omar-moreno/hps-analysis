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
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TF1.h>

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <RootFileReader.h>

class ComparePlots : public RootFileReader { 

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
        void overlayPlots();

        /**
         *
         */
        void setStyle(const std::string style) { this->style = style; }; 

    private: 

        /**
         *
         */
        void applyBasic1DStyle();

        /**
         *
         */
        void applyBasic2DStyle();

        /**
         *
         */
        void applyMCStyle();

        std::string style;
};

#endif // __COMPARE_PLOTS_H__

