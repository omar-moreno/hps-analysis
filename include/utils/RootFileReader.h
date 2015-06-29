/**
 *  @file RootFileReader.h
 *  @brief Reader used to parse ROOT files. 
 *  @author Omar Moreno <omoreno1@ucsc.edu>
 *  @date June 22, 2015
 */

#ifndef __ROOT_FILE_READER_H__
#define __ROOT_FILE_READER_H__

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <iostream>
#include <map>
#include <vector>
#include <list>

//------------//
//--- ROOT ---//
//------------//
#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TGraph.h>

class RootFileReader { 

    public:

        /**
         * Default Constructor
         */
        RootFileReader();

        /**
         * Destructor
         */
        ~RootFileReader();

        /**
         * Open a ROOT file, parse it and load all histograms to their 
         * corresponding maps.
         *
         * @param root_file ROOT file to parse
         */
        void parseFile(TFile* root_file);

        /**
         *
         */
        void parseFile(TFile* root_file, std::string histogram_substring);

        /**
         *
         */ 
        void parseFiles(std::list<TFile*> root_files);

        /**
         *
         */
        void parseFiles(std::list<TFile*> root_files, std::string histogram_substring);


        /**
         *
         */
        void setHistogramName(std::string histogram_substring) { this->histogram_substring = histogram_substring; };

        /**
         *
         *
         */
        std::vector<TH1*> getMatching1DHistograms(std::string histogram_name);

        /**
         *
         */
        std::vector<TH1*> getMatching2DHistograms(std::string histogram_name);

        /**
         *
         */
        std::vector<TGraph*> getMatchingGraphs(std::string histogram_name);

    protected: 

        /** Map containing 1D histograms */
        std::map <std::string, std::vector<TH1*> > histogram1D_map;
        
        /** Map containing 2D histograms */
        std::map <std::string, std::vector<TH1*> > histogram2D_map;

        /** Map containing objects of type TGraph */
        std::map <std::string, std::vector<TGraph*> > graph_map;

    private:

        /**
         * Parse and load all histograms to their corresponding maps.
         *
         * @param keys List of keys retrieved from a ROOT file
         */
        void parseFile(TList* keys);

        /** Name (or part) of the histogram of interest in stored in the file */
        std::string histogram_substring;

}; // RootFileReader

#endif // RootFileReader
