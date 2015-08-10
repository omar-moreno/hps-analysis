/**
 * @file GblTrackAnalysis.h
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date August 4, 2015
 *
 */

#ifndef __GBL_TRACK_ANALYSIS_H__
#define __GBL_TRACK_ANALYSIS_H__

#include <iostream>

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <HpsAnalysis.h>
#include <Plotter.h>
#include <TrackExtrapolator.h>

//---------------//
//--- HPS DST ---//
//---------------//
#include <GblTrack.h>
#include <SvtTrack.h>

class GblTrackAnalysis : public HpsAnalysis { 

    public: 
       
        /**
         * Constructor
         */
        GblTrackAnalysis();

        /**
         *  Destructor
         */
        ~GblTrackAnalysis();

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

        Plotter* plotter;

        // Name of the class
		std::string class_name;


}; // GblTrackAnalysis



#endif // __GBL_TRACK_ANALYSIS_H__
