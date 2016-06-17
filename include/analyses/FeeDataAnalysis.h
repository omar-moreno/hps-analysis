/**
 *
 * @file FeeDataAnalysis.h
 * @brief Analysis used to study Coulomb scattered electrons using HPS data.
 * @author Omar Moreno <omoreno@slac.stanford.edu>
 *         SLAC
 * @date June 14, 2016
 * 
 */

#ifndef __FEE_DATA_ANALYSIS_H__
#define __FEE_DATA_ANALYSIS_H__

//------------------//
//   HPS Analysis   //
//------------------//
#include <FeeAnalysis.h>

class FeeDataAnalysis : public FeeAnalysis { 
    
    /**
     * Process an HPS event i.e. {@link HpsEvent} object. 
     *
     * @param event {@link HpsEvent} object to process.
     */
    void processEvent(HpsEvent* event);

    /** @return A string representation of this analysis. */
    std::string toString(); 

};

#endif
