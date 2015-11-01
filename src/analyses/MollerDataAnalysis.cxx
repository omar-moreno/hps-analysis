/**
 *
 * @file MollerDataAnalysis.cxx
 * @brief Analysis used to select Moller events  in the Engineering Run 2015 
 *        data.
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date October 30, 2015
 * 
 */

#include <MollerDataAnalysis.h>

void MollerDataAnalysis::processEvent(HpsEvent* event) { 

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
    MollerAnalysis::processEvent(event); 
}

std::string MollerDataAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
