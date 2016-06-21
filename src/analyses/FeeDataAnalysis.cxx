/**
 *
 * @file FeeDataAnalysis.cxx
 * @brief Analysis used to study Coulomb scattered electrons using HPS data.
 * @author Omar Moreno <omoreno@slac.stanford.edu>
 *         SLAC
 * @date June 14, 2016
 * 
 */

#include <FeeDataAnalysis.h>

FeeDataAnalysis::FeeDataAnalysis() 
    : class_name("FeeDataAnalysis") { 
}

void FeeDataAnalysis::processEvent(HpsEvent* event) { 

    // Only look at pairs1 triggers
    if (!event->isSingle1Trigger()) return;

    // Only look at events with the SVT bias ON
    if (!event->isSvtBiasOn()) return; 

    // Only look at events where the SVT is closed
    if (!event->isSvtClosed()) return;

    // Skip events that had burst mode noise
    if (event->hasSvtBurstModeNoise()) return; 

    // Skip events that had SVT header errors
    if (event->hasSvtEventHeaderErrors()) return;

    // 
    FeeAnalysis::processEvent(event); 
}

std::string FeeDataAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
