/**
 * @file SharedHitAnalysis.h
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date June 15, 2015
 */

#ifndef __SHARED_HIT_ANALYSIS_H__
#define __SHARED_HIT_ANALYSIS_H__

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <HpsAnalysis.h>
#include <Plotter.h>

//---------------//
//--- HPS DST ---//
//---------------//
#include <SvtTrack.h>
#include <SvtHit.h>

class SharedHitAnalysis : public HpsAnalysis { 

    public: 

        /**
         * Constructor
         */
        SharedHitAnalysis();

        /**
         * Destructor
         */
        ~SharedHitAnalysis();

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

        SvtTrack* first_track;
        SvtTrack* second_track;

        Plotter* track_plotter;

        // Name of the class
		std::string class_name;

}; // SharedHitAnalysis

#endif // __SHARED_HIT_ANALYSIS_H__
