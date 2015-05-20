/**
 * @file TrackAnalysis.h
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date May 14, 2015
 *
 */

#ifndef __TRACK_ANALYSIS_H__
#define __TRACK_ANALYSIS_H__

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <unordered_map>
#include <vector>

//------------//
//--- ROOT ---//
//------------//
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <HpsAnalysis.h>

//---------------//
//--- HPS DST ---//
//---------------//
#include <SvtTrack.h>
#include <SvtHit.h>

class TrackAnalysis : public HpsAnalysis { 

    public: 
       
        /**
         * Constructor
         */
        TrackAnalysis();

        /**
         *  Destructor
         */
        ~TrackAnalysis();

        /**
         *  Method used to initialize an HPS analysis.
         */
        void initialize();

        /**
         *  Method containing the code used to process an HpsEvent.
         *
         *  @param event : HpsEvent that will be processed
         */
        void processEvent(HpsEvent* event);

        /**
         *  Method used to finalize an HPS analysis.
         */
        void finalize();

        /**
         *  Method used to initialize any histograms used by the analysis.
         */
        // TODO:  This should use a histogram factory instead
        void bookHistograms();

        /**
         *  Provide a string representation of this analysis.
         *
         *  @return String representation of this analysis.
         */
        std::string toString(); 
   
    private:

        /**
         *
         */
        double getMagnitude(std::vector<double> v);

        SvtTrack* track;

        TCanvas* canvas; 
        
        // ROOT output_file
        TFile* output_file; 

        std::unordered_map<std::string, TH1F*> track_parameter_plots;
        std::unordered_map<std::string, TH1F*> track_plots;
        std::unordered_map<std::string, TH1F*> electron_track_plots;
        std::unordered_map<std::string, TH1F*> positron_track_plots;

        std::unordered_map<std::string, TH2F*> hit_position_plots;
        std::unordered_map<std::string, TH2F*> track_plots_2d; 
        std::unordered_map<std::string, TH2F*> track_epem_2d_plots; 
        std::unordered_map<std::string, TH2F*> track_emem_2d_plots; 

        // Name of the class
		std::string class_name;


}; // TrackAnalysis

#endif // __TRACK_ANALYSIS_H__
