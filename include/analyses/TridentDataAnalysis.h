/**
 *
 * @file TridentDataAnalysis.h
 * @brief Analysis used to study Tridents in the Engineering Run 2015 data.
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date October 14, 2015
 *
 */

#ifndef __TRIDENT_DATA_ANALYSIS_H__
#define __TRIDENT_DATA_ANALYSIS_H__

//------------------//
//   HPS Analysis   //
//------------------//
#include <TridentAnalysis.h>

class TridentDataAnalysis : public TridentAnalysis { 

    /**
     * Process an HPS event i.e. {@link HpsEvent} object and extract Trident
     * candidates from the data.
     *
     * @param event {@link HpsEvent} object to process.
     */
    void processEvent(HpsEvent* event);

    /** @return A string representation of this analysis. */
    std::string toString(); 

}; // TridentDataAnalysis

#endif // __TRIDENT_DATA_ANALYSIS_H__
