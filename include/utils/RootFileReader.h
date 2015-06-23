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

//------------//
//--- ROOT ---//
//------------//
#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TH1.h>

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
         *
         */
        std::vector<TH1*> getMatching1DHistograms(std::string histogram_name);

    private:

        /**
         * Parse and load all histograms to their corresponding maps.
         *
         * @param keys List of keys retrieved from a ROOT file
         */
        void parseFile(TList* keys);

        /** Map containing 1D histograms */
        std::map <std::string, std::vector<TH1*> > histogram1D_map;
        
        /** Map containing 2D histograms */
        std::map <std::string, std::vector<TH1*> > histogram2D_map;

}; // RootFileReader

#endif // RootFileReader
