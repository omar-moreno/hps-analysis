/**
 *
 * @file V0Analysis.h
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date September 29, 2015
 * @brief Analysis that uses V0 particles to look at Tridents.
 *
 */

#ifndef __V0_ANALYSIS_H__
#define __V0_ANALYSIS_H__

//------------------//
//   HPS Analysis   //
//------------------//
#include <EcalUtils.h>
#include <HpsAnalysis.h>
#include <FlatTupleMaker.h>
#include <TrackClusterMatcher.h>

class V0Analysis : public HpsAnalysis { 

    public: 

        /** Constructor */
        V0Analysis(); 

        /** Destructor */
        ~V0Analysis();
        
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

        /** A set of Ecal utilities */
        EcalUtils* ecal_utils;

        /** */
        TrackClusterMatcher* matcher; 

        /** */
        FlatTupleMaker* tuple; 

        /** Name of the class */
        std::string class_name; 
};



#endif // __V0_ANALYSIS_H__
