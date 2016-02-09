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

//------------------//
//   HPS Analysis   //
//------------------//
#include <EcalUtils.h>
#include <HpsAnalysis.h>
#include <FlatTupleMaker.h>
#include <AnalysisUtils.h>
#include <TrackClusterMatcher.h>

//-------------//
//   HPS DST   //
//-------------//
#include <SvtTrack.h>

class MollerAnalysis : public HpsAnalysis { 

    public: 

        /** Constructor */
        MollerAnalysis();

        /** Destructor */
        ~MollerAnalysis();
        
        /** Initialize an HPS analysis. */
        void initialize();

        /**
         * Process an HpsEvent.
         *
         * @param event HpsEvent that will be processed.
         */
        void processEvent(HpsEvent* event);

        /** Finalize an HPS analysis. */
        void finalize();

        /** Initialize histograms used in this analysis. */
        void bookHistograms();

        /** @return A string representation of this analysis. */
        std::string toString();
   
    protected: 

        /** Name of the class */
        std::string class_name; 

    private:

        /** A set of Ecal utilities */
        EcalUtils* ecal_utils;

        /** */
        TrackClusterMatcher* matcher; 

        /** */
        FlatTupleMaker* tuple; 
 
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

        /** Total number of events with a "good" cluster pair. */
        int cluster_pair_counter;  
};

#endif // __MOLLER_ANALYSIS_H__
