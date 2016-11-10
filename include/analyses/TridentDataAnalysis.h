/**
 *
 * @file TridentDataAnalysis.h
 * @brief Analysis used to study Tridents in the Engineering Run 2015 data.
 * @author Omar Moreno <omoreno@slac.stanford.edu>
 *         SLAC National Accelerator Laboratory
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

    public: 

        /**
         * Process an HPS event i.e. {@link HpsEvent} object and extract Trident
         * candidates from the data.
         *
         * @param event {@link HpsEvent} object to process.
         */
        void processEvent(HpsEvent* event);

        /** Finalize an HPS analysis. */
        void finalize(); 

        /** @return A string representation of this analysis. */
        std::string toString(); 

    private: 
    
        /** Trigger count */
        double trigger_count{};

        /** Good SVT event count */
        double svt_quality_count{};


}; // TridentDataAnalysis

#endif // __TRIDENT_DATA_ANALYSIS_H__
