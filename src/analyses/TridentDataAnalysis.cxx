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

    // Only look at pair1 triggers.
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
    ++svt_quality_count;

    // Use the base class to process the event. 
    TridentAnalysis::processEvent(event); 
}

void TridentDataAnalysis::finalize() { 
    std::cout << std::fixed;
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    std::cout << "%   Data Only  " << std::endl;
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    std::cout << "% Trigger count: " << trigger_count << std::endl;
    std::cout << "% SVT quality count: " << svt_quality_count << std::endl;
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;

    TridentAnalysis::finalize();
}

std::string TridentDataAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
