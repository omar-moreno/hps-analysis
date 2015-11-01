/**
 *
 * @file MollerDataAnalysis.h
 * @brief Analysis used to select Moller events  in the Engineering Run 2015 
 *        data.
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date October 30, 2015
 * 
 */

#ifndef __MOLLER_DATA_ANALYSIS_H__
#define __MOLLER_DATA_ANALYSIS_H__

//------------------//
//   HPS Analysis   //
//------------------//
#include <MollerAnalysis.h>

class MollerDataAnalysis : public MollerDataAnalysis { 

    /**
     * Process an HPS event i.e. {@link HpsEvent} object and extract Moller
     * candidates from the data.
     *
     * @param event {@link HpsEvent} object to process.
     */
    void processEvent(HpsEvent* event);

    /** @return A string representation of this analysis. */
    std::string toString(); 


}; // MollerDataAnalysis 

#endif
