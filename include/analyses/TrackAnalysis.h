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

//------------//
//--- ROOT ---//
//------------//
#include <TCanvas.h>
#include <TH1F.h>

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <HpsAnalysis.h>

//---------------//
//--- HPS DST ---//
//---------------//
#include <SvtTrack.h>

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

        SvtTrack* track;

        TCanvas* canvas; 

        std::unordered_map<std::string, TH1F*> track_parameter_plots;

}; // TrackAnalysis

#endif // __TRACK_ANALYSIS_H__
