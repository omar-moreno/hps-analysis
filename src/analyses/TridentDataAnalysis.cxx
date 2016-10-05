/**
 *
 * @file TridentDataAnalysis.h
 * @brief Analysis used to study Tridents in the Engineering Run 2015 data.
 * @author Omar Moreno <omoreno@slac.stanford.edu>
 *         SLAC National Accelerator Laboratory
 * @date October 14, 2015
 *
 */

#include <TridentDataAnalysis.h>

void TridentDataAnalysis::processEvent(HpsEvent* event) { 

    // Only look at pairs1 triggers.
    if (!event->isPair1Trigger()) return;
    ++trigger_count;

    // Only look at events with the SVT bias ON.
    if (!event->isSvtBiasOn()) return; 

    // Only look at events where the SVT is closed.
    if (!event->isSvtClosed()) return;

    // Skip events that had burst mode noise.
    if (event->hasSvtBurstModeNoise()) return; 

    // Skip events that had SVT header errors.
    if (event->hasSvtEventHeaderErrors()) return;
    ++svt_event_count;

    // Use the base class to process the event. 
    TridentAnalysis::processEvent(event); 
}

std::string TridentDataAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
