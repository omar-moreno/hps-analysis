
#include <V0DataAnalysis.h>

void V0DataAnalysis::processEvent(HpsEvent* event) { 
    
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
    V0Analysis::processEvent(event); 

}

std::string V0DataAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
