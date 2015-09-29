/**
 * @file MollerAnalysis.h
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date July 10, 2015
 */

#ifndef __MOLLER_ANALYSIS_H__
#define __MOLLER_ANALYSIS_H__

#include <cmath>

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <HpsAnalysis.h>
#include <Plotter.h>
#include <AnalysisUtils.h>
#include <TrackExtrapolator.h>
#include <TrackClusterMatcher.h>
#include <EcalUtils.h>

//---------------//
//--- HPS DST ---//
//---------------//
#include <SvtTrack.h>

class MollerAnalysis : public HpsAnalysis { 

    public: 

        /**
         * Constructor
         */
        MollerAnalysis();

        /**
         * Destructor
         */
        ~MollerAnalysis();

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

        // TODO: Move this class to a utility class
        std::vector<EcalCluster*> getClusterPair(HpsEvent* event);

        bool isMatch(EcalCluster* cluster, SvtTrack* track);
        
        Plotter* plotter;
        
        TrackClusterMatcher* matcher;

        EcalUtils* ecal_utils; 

        // Name of the class
        std::string class_name;

        //-- Event counters --//
        //--------------------//
        
        /** Total number of events */
        int event_counter;

        /** Total number of singles1 trigger events where the SVT bias was on */
        int bias_on_counter; 

        /** Total number of singles1 triggers */
        int single1_trigger_counter;

        /** 
         * Total number of singles 1 triggers where the SVT bias was on and it
         * was in closed position.
         */
        int svt_closed_position_counter; 
};

#endif // __MOLLER_ANALYSIS_H__
