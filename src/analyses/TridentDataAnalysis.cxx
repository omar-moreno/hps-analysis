/**
 *
 * @file TridentDataAnalysis.cxx
 * @brief Analysis used to study Tridents in the Engineering Run 2015 data.
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date October 14, 2015
 *
 */

#include <TridentDataAnalysis.h>

void TridentDataAnalysis::processEvent(HpsEvent* event) { 

    // Only look at pairs1 triggers
    if (!event->isPair1Trigger()) return;

    // Only look at events with the SVT bias ON
    if (!event->isSvtBiasOn()) return; 

    // Only look at events where the SVT is closed
    if (!event->isSvtClosed()) return;

    // Skip events that had burst mode noise
    if (event->hasSvtBurstModeNoise()) return; 

    // Skip events that had SVT header errors
    if (event->hasSvtEventHeaderErrors()) return;

    // 
    TridentAnalysis::processEvent(event); 
}

std::string TridentDataAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
