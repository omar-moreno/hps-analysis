/**
 * @file TagProbeDataAnalysis.h
 * @brief Analysis used to calculate tracking efficiency using Moller events
 *        in the Engineering Run 2015 data.
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date February 16, 2015
 */

#ifndef __TAG_PROBE_DATA_ANALYSIS_H__
#define __TAG_PROBE_DATA_ANALYSIS_H__

#include <TagProbeAnalysis.h>

class TagProbeDataAnalysis : public TagProbeAnalysis {

    public: 
    
        /**
         * Process an HPS event i.e. {@link HpsEvent} object and extract Moller
         * candidates from the data.
         *
         * @param event {@link HpsEvent} object to process.
         */
        void processEvent(HpsEvent* event);

        /** @return A string representation of this analysis. */
        std::string toString() { return "TagProbeDataAnalysis"; }; 

}; // TagProbeDataAnalysis 

#endif // __TAG_PROBE_DATA_ANALYSIS_H__
